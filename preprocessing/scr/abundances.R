# Abundance and compositional analysis #########################################

#panImm,Tcell,validation	comp	abundances	abundances
#panImm,Tcell,validation	ords	abundances	abundances

#panImm	da_res	preprocess	differential_abundance
#PC_scores	da_res	preprocess	differential_abundance
#Tcell	da_res	preprocess	differential_abundance
#Tcell_obs	abundances	preprocess	differential_abundance
#Valid_obs	abundances	preprocess	differential_abundance


# Functions ####################################################################

# Compositional data analysis (CODA) framework
comp.prep = function(ab.in,celltype){
  
  ab = ab.in %>% as_tibble() %>%
    pivot_wider(names_from=celltype,values_from='n',values_fill=0) %>%
    pivot_longer(!sample_id,names_to='celltype',values_to='n') %>%
    group_by(sample_id) %>% mutate(prop = n / sum(n)) %>% ungroup()
  
  #CODA
  ab %>% select(-prop) %>%
    pivot_wider(names_from=celltype,values_from='n',values_fill=0) -> ab.mat
  #cmultRepl handle zeroes, apply CLR transformation
  as_tibble(zCompositions::cmultRepl(select(ab.mat,-sample_id), 
                                     output = "p-counts", z.warning=0.99)) %>%
    mutate(sample_id = ab.mat$sample_id) %>%
    pivot_longer(!sample_id,names_to='celltype',values_to='pcount') %>%
    group_by(sample_id) %>%
    mutate(pcount.prop = pcount / sum(pcount)) %>% #pseudocount proportions
    mutate(clr = as.vector(compositions::clr(pcount))) %>%
    ungroup() -> ab.clr
  #comp data
  comp = ab %>% left_join(ab.clr %>% select(sample_id,celltype,clr)) %>% 
    left_join(
      #summary data
      ab %>% group_by(sample_id) %>% summarise(sampleSize=sum(n)) %>% 
        left_join(ab.clr %>% group_by(sample_id) %>% 
                    summarise(geom.mean = exp(mean(log(pcount.prop))))
      )
  )
  comp
}

# Ordination
comp.ord = function(comp.){
  
  # For celltypes detected in min 5% of samples
  celltype_detected = comp. %>% group_by(celltype) %>%
    summarise(pctGrt0=sum(n>0)/n()*100, mean=mean(clr), med=median(clr), sd=sd(clr)) %>%
    arrange(pctGrt0) %>% filter(celltype>5) %>% pull(celltype)
  
  mat = comp. %>%
    filter(celltype %in% celltype_detected) %>%
    select(sample_id,celltype,clr) %>%
    pivot_wider(names_from=celltype, values_from=clr) %>%
    data.frame() %>% column_to_rownames('sample_id') %>%
    scale()
  
  # Ordination
  ord <- vegan::rda(mat)
  
  list(mat=mat, ord=ord)
}



# Import #######################################################################

panImm = read.csv('data/panImm_obs.csv')
Tcell = read.csv('data/Tcell_obs.csv')
valid = read.csv('../integration/data/validation_Tcell-obs.csv') %>% 
  left_join(read.csv('../integration/data/validation_Tcell.csv')) %>% 
  filter(!pheno_pred %in% c('Unknown','Myeloid','CD3_neg','Doublet')) %>% 
  filter(pheno_uncertainty < 0.2)
valid_md = as_tibble(read.csv('../integration/data/validation-metadata.csv'))

# Compositional data ###########################################################

comp = list()

md_cols = c('sample_id','study_id','donor_id','tissue','chem','sort')

comp[['panImm']] = panImm %>% 
  group_by(sample_id,lineage) %>% tally() %>% 
  comp.prep(.,'lineage') %>% 
  left_join(panImm %>% select(all_of(md_cols)) %>% distinct()) %>%
  left_join(read.csv('data/import/metadata-donor.csv') %>% distinct(),
            relationship = "many-to-many")
  
comp[['Tcell']] = Tcell %>% 
  group_by(sample_id,pheno) %>% tally() %>% 
  comp.prep(.,'pheno') %>% 
  left_join(Tcell %>% select(all_of(md_cols)) %>% distinct()) %>%
  left_join(read.csv('data/import/metadata-donor.csv') %>% distinct(),
            relationship = "many-to-many")
  
comp[['valid']] = valid %>% 
  group_by(sample_id,pheno_pred) %>% tally() %>% 
  comp.prep(.,'pheno_pred') %>% 
  left_join(Tcell %>% select(all_of(md_cols)) %>% distinct()) %>%
  left_join(valid_md)

# Ordination (PCA) #############################################################

ords=list()

# Remove selected donors
donors.longit = c('Foster_2024.MM4',"Foster_2024.MM3",'Foster_2024.SMM13',
                  'Liu_2021.37692_SMM','Liu_2021.47491_SMM',
                  'Liu_2021.58408_SMM','Oetjen_2018.C_T2','Oetjen_2018.S_T2')
donor.Tex_hi = 'Zheng_2021.P20190122'
donors.remove = c(donors.longit,donor.Tex_hi)

ords[['Tcell']] = comp$Tcell %>% 
  filter(sampleSize>100, !donor_id %in% donors.remove, tissue=='BM') %>% 
  comp.ord() %>% .$ord

ords[['Tcell']]$PCA = ords[['Tcell']]$CA$u[,1:5] %>% data.frame() %>% rownames_to_column('sample_id') %>% as_tibble()

# Validation data PCA re-calculation ###########################################

comp.dat = comp$valid %>% filter(!celltype %in% c('CD8.Tex','CD8.Trm'))

# Original Tcell ord matrix
mat.orig.dim = as.matrix(ords$Tcell[["Ybar"]]) %>% dim()
mat.orig = as.matrix(ords$Tcell[["Ybar"]])[1:mat.orig.dim[1],1:mat.orig.dim[2]]

# Remove missing pheno
mat.orig=mat.orig[,colnames(mat.orig) %in% comp.dat$celltype]
comp.dat = comp.dat[comp.dat$celltype %in% colnames(mat.orig),]

# Validation data matrix, scaled to center/range of original matrix  
mat_valid = comp.dat %>%
  select(sample_id,celltype,clr) %>%
  pivot_wider(names_from=celltype, values_from=clr) %>%
  data.frame() %>% column_to_rownames('sample_id') %>% 
  scale(., center = colMeans(mat.orig), scale = matrixStats::colSds(mat.orig))
  
# Re-calculatr PCA, matrix multiplication 
pca_valid = as.matrix(mat_valid) %*% ords$Tcell$CA$v[colnames(mat_valid),] %>%
  data.frame() %>% .[,1:2] %>%
  rownames_to_column('sample_id') %>% as_tibble()

ords[['valid']] = pca_valid

# Save outputs #################################################################

save(comp,file="data/comp.RData")
save(ords,file="data/ords.RData")
