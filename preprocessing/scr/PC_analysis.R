# Plasma cell analysis  ########################################################

library(tidyverse)
library(zellkonverter)
library(scran)
library(UCell)
library(msigdbr)
library(fgsea)

# Plasma cell pathway scores ###################################################

sce = readH5AD('data/import/PC.h5ad')
mPC_id = read.csv('data/mPC_id.csv') %>% filter(mPC=='True')
pathways = read.csv('../resources/pathways.csv')
gs = list()
for ( g in unique(pathways$pathway)[1:40]){gs[[g]]=pathways[pathways$pathway==g,]$gene}

# UCell scoring
sce.mat = sce@assays@data@listData[["X"]]
rownames(sce.mat) = rownames(sce)
ucell.scores <- UCell::ScoreSignatures_UCell(sce.mat, features=gs)
ucell.scores = data.frame(ucell.scores)
colnames(ucell.scores) = str_remove(colnames(ucell.scores), '_UCell')
scr = SingleCellExperiment(assay=list(scores=t(ucell.scores)), colData=colData(sce))
save(scr, file='data/tumour_modules_scr.RData')

# Plasma cell pathway scores, high/low quantification ##########################

load('data/tumour_modules_scr.RData')
pathways = read.csv('../resources/pathways.csv')[,1:2] %>% distinct() %>% .[1:40,]
rowData(scr) = DataFrame(data.frame(
  p=rownames(scr),pathway=pathways$pathway,pathway_neat=pathways$pathway_neat))
mPC_id = read.csv('data/mPC_id.csv') %>% filter(mPC=='True')

# Patients with sufficient tumour cells
full = as.data.frame(scater::makePerCellDF(
  scr,rownames(scr),exprs_values='scores')) %>%
  filter(colnames(scr) %in% mPC_id$obs_names)
full = full %>% filter(donor_id %in% names(table(full$donor_id) %>% .[.>30]))

# For each pathway, calculate the fraction of cells 1sd > median
# (highly expressing said pathway)
res.list = list()
for ( mod in rownames(scr) ){
  # Expression
  full$pl = full[,mod]
  # Pathway thr
  thr = median(full$pl)+sd(full$pl)
  # Quantify % above thr
  res.list[[mod]] = full %>% mutate(is = pl > thr) %>%
    group_by(donor_id,is) %>% tally() %>% ungroup() %>%
    group_by(donor_id) %>% mutate(pct=n/sum(n)*100) %>% ungroup() %>%
    filter(is) %>% select(donor_id,pct)
}

bind_rows(res.list,.id='p') %>% as_tibble() %>%
  left_join(data.frame(rowData(scr))) %>%
  write.csv('data/tumour_modules_pct.csv',row.names=F)

# Tumour/non-tumour DEG ########################################################

sce = readH5AD('data/import/PC.h5ad')
mPC_id = read.csv('data/mPC_id.csv')
Foster_2024_batch_id = read.csv('../resources/Foster_2024_batch_id.csv')

# mPC calls
sce$mPC = colnames(sce) %in% mPC_id[mPC_id$mPC=='True',]$obs_names

# Subset to donors/batches
# Donors w mPC
mPC_donor = unique(mPC_id[mPC_id$mPC=='True',]$donor)
# Batches where donors had mPC identified
# (Pre-filtered for sufficient mPC)
mPC_batch = as_tibble(colData(sce)) %>%
  filter(donor_id %in% mPC_donor, !mPC) %>% pull(study_id) %>% unique()
# Control donors
donor_Non = read.csv('data/import/metadata-donor.csv') %>%
  .[.$cohort=='Non',] %>% pull(donor_id)
# Final subset (used Foster_2024 batch1 and batch2)
Foster_2024_use = Foster_2024_batch_id[Foster_2024_batch_id$batch_id %in%
                                         c('batch1','batch2'),]$donor_id
index.use = as_tibble(colData(sce)) %>% mutate(index=colnames(sce)) %>%
  filter(donor_id %in% Foster_2024_use) %>%
  filter( donor_id %in% mPC_donor | donor_id %in% donor_Non) %>% pull(index)

# Subset to non-proliferating non-ISG+ (using pathway scores)
colData(sce) = as_tibble(colData(sce)) %>%
  mutate(index=colnames(sce)) %>%
  left_join(
    as_tibble(data.frame(index=colnames(scr),
      Cycle=as.vector(assay(scr['MP2..Cell.Cycle...G1.S',])),
      ISG=as.vector(assay(scr['MP17.Interferon.MHC.II..I.',]))))
  ) %>% DataFrame()
#as_tibble(colData(sce)) %>% ggplot(aes(Cycle))+geom_histogram(bins=100)+geom_vline(xintercept=0.075)
#as_tibble(colData(sce)) %>% ggplot(aes(ISG))+geom_histogram(bins=100)+geom_vline(xintercept=0.4)

# Subset cells
cells.use = as_tibble(colData(sce)) %>%
  filter(index %in% index.use, Cycle<0.075, ISG< 0.4) %>% pull(index)
sce = sce[,sce$index %in% cells.use]

# donor_id/batch_id
colData(sce) = as_tibble(colData(sce)) %>%
  left_join(Foster_2024_batch_id) %>% DataFrame()
donor_batch = as_tibble(colData(sce)) %>%select(donor_id,batch_id) %>% distinct()
donor_batch = setNames(donor_batch$batch_id,donor_batch$donor_id)

# Iterative testing over donors
donor.res = list()
for ( donor.i in intersect(mPC_donor,unique(sce$donor_id)) ){
  # Subset to this donor OR this donor's batch non-clonal PC
  donor.bool = sce$donor_id==donor.i | ( sce$batch_id==donor_batch[[donor.i]] & !sce$mPC )
  sce.i = sce[,donor.bool]
  # Test between mPC (for this donor only) and non-mPC (for entire batch)
  donor.res[[donor.i]] = scran::pairwiseTTests(sce.i,
    groups=sce.i$mPC, assay.type='X')$statistics[[1]] %>%
      data.frame() %>% rownames_to_column('gene') %>% as_tibble() %>%
      filter(FDR<0.1) %>% select(gene,FDR,logFC)
}

# Remove certain gene groups, remove longitudinal donors
ig_genes = read.csv('../resources/ig_genes.csv')$gene
ribo_genes = read.csv('../resources/ribo_genes.csv')$gene
mt_genes = as_tibble(rowData(sce)) %>% filter(mt) %>% pull(gene)
donors.longit = c('Foster_2024.MM4',"Foster_2024.MM3",'Foster_2024.SMM13')

tumour_DE = bind_rows(donor.res, .id='donor_id') %>%
  # Negative FC = UP tumour
  mutate(logFC=logFC*-1, de=ifelse(logFC>0,'tum','non')) %>%
  # Remove logit donors, need md
  filter(!donor_id %in% donors.longit) %>%
  # Remove noise-associated gene groups
  filter(!gene %in% c(ig_genes,mt_genes,ribo_genes))

tumour_DE %>% write.csv('data/tumour_DE.csv')

# Tumour/non-tumour DEG pathway analysis #######################################

# GSEA gene sets
gs=list()
sets <- c("CP:BIOCARTA","CP:KEGG","CP:REACTOME")
msig <- subset(msigdbr(species = "human", category = "C2"), gs_subcat %in% sets)
for (geneset in unique(msig$gs_name)){gs[[geneset]]<-subset(msig, gs_name==geneset)$gene_symbol}

# Hallmarks
msig = msigdbr(species = "human", category = "H")
for (geneset in unique(msig$gs_name)){gs[[geneset]]<-subset(msig, gs_name==geneset)$gene_symbol}

# Gavish et al. genesets
for (geneset in unique(pathways$pathway)[1:40]){gs[[geneset]]=pathways[pathways$pathway==geneset,]$gene}

# fgsea over each individual's tumour cells
fgseas = list()
for ( i in unique(tumour_DE$donor_id)){
  ranks = tumour_DE %>% filter(donor_id==i) %>% select(gene,logFC) %>% deframe()
  fgseas[[i]] = as_tibble(fgsea(gs, ranks, maxSize=500)) %>% filter(padj<0.1)
}
fgseas = bind_rows(fgseas, .id='donor_id') %>% mutate(de=ifelse(NES>0,'tum','non'))
save(fgseas,file='data/tumour_DE_fgsea.RData')
























