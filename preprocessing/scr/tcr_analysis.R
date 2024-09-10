# TCR analysis #################################################################

library(tidyverse)
# Simpson's diversity calculation code
simpson.sample=function(x){
  n=sum(x,na.rm=TRUE)
  pc=sum((x-1)*x/(n*(n-1)))
  sp_d=1/(pc)
  sp_d
}
source('../resources/aes.R')

obs = read.csv('data/Tcell_obs.csv')
tcr = read.csv('data/import/Tcell-scTCR.csv')
donor_md = read.csv('data/import/metadata-donor.csv')

tcr_analysis = list()

# CD4/CD8 subset identification ################################################

# Count N CD4+ and CD8+ cells per clone, threshold to two subsets and DN
# (double negative, unknown)

tcr.subset = tcr %>% group_by(clone_id) %>%
  summarise(n=n(),CD8_bool=sum(CD8_Ncells),CD4_bool=sum(CD4_Ncells)) %>%
  mutate(CD8_prop=CD8_bool/n, CD8_thr=CD8_prop>0.25,CD8_q=paste0('Quantile ',ntile(CD8_prop,4)))

tcr.subset = tcr.subset %>% mutate(type = case_when(
  n==1 & CD8_bool==0 & CD4_bool==0 ~ 'Unknown',
  CD4_bool<CD8_bool ~ 'CD8',
  CD4_bool>CD8_bool ~ 'CD4',
  CD4_bool==CD8_bool ~ 'CD4',
))

tcr_analysis[['type']] = tcr.subset %>% select(clone_id,type)


# Clonality ####################################################################

# Calculate clonality across all cells and among CD8+ memory T cells (Tm)
CD8.Tm = c("CD8.Tem.IL7R","CD8.Tem.KLRG1","CD8.TemActive",
           "CD8.TEMRA","CD8.TemTerm","CD8.Tex","CD8.Trm")

div.list = list()

# All T cells
div.list[['all']] = obs %>% filter(clone_id!='') %>% 
  group_by(donor_id,clone_id) %>% tally() %>% ungroup() %>%
  group_by(donor_id) %>% summarise(
    counts=sum(n),unique=length(unique(clone_id)),div=simpson.sample(n)
  ) %>% ungroup() %>%
  filter(counts>100)

# CD8+ clones, CD8+Tm cluster
div.list[['CD8.Tm']] = obs %>% filter(clone_id!='') %>% 
  filter(clone_id %in% tcr.subset[tcr.subset$type=='CD8',]$clone_id) %>%
  filter(pheno %in% CD8.Tm) %>%
  group_by(donor_id,clone_id) %>% tally() %>% ungroup() %>%
  group_by(donor_id) %>% summarise(
    counts=sum(n),unique=length(unique(clone_id)),div=simpson.sample(n)
  ) %>% ungroup() %>%
  filter(counts>100)

# Results
tcr_analysis[['clonality']] = bind_rows(div.list, .id='group') %>%
  # for ease of plotting
  mutate(clon=log10(1/div)) %>%
  left_join(donor_md %>% select(donor_id,age,cohort))

# TCR clustering ###############################################################

tcr_analysis[['clus']] = list()

# Import tcrdist3 results
tcr_analysis$clus[['clusters']] = read_csv('data/tcr_clustering.csv')

# metaclone ID
tcr_analysis$clus[['clus_id']] = tcr_analysis$clus[['clusters']] %>%
  select(subject,chain,cdr3s) %>% distinct() %>%
  group_by(subject) %>% mutate(clus_id = paste0(subject,'-',chain,'-',row_number())) %>%
  ungroup()

# clus id to clone_id
tcr_analysis$clus[['clon_clus']] = tcr_analysis$clus[['clusters']] %>% 
  left_join(tcr_analysis$clus$clus_id) %>%
  mutate(v_gene = str_replace(v_gene, '\\*01', ''), j_gene = str_replace(j_gene, '\\*01', '')) %>%
  rename(donor_id=subject) %>%
  left_join(tcr %>% select(donor_id,clone_id,CDR3aa,v_gene,j_gene) %>% distinct()) %>%
  select(clus_id,clone_id) %>% distinct()

# N clones per cluster
tcr_analysis$clus[['Nclone']] = tcr_analysis$clus$clon_clus %>% group_by(clus_id) %>% tally()
tcr_analysis$clus[['Nclone']]$exp = tcr_analysis$clus[['Nclone']]$n > 1

# N/pct cells in pheno per cluster
tcr_analysis$clus[['pheno']] = obs %>%
  filter(clone_id %in% tcr_analysis$clus$clon_clus$clone_id) %>%
  group_by(donor_id,clone_id,pheno) %>% tally() %>% ungroup() %>%
  group_by(donor_id,clone_id) %>% mutate(pct=n/sum(n)*100) %>%
  left_join(tcr_analysis$clus$clon_clus)

# % repertoire clustered
clon_clus.sub = tcr_analysis$clus$clon_clus %>%
  left_join(tcr_analysis$clus$clus_id %>% select(subject,chain,clus_id)) %>%
  distinct() %>% rename(donor_id=subject)

clustered = obs %>% select(donor_id,clone_id,tissue) %>%
  filter(clone_id %in% (tcr %>% filter(clone_size>1) %>% pull(clone_id))) %>%
  filter(tissue=='BM') %>% select(-tissue) %>%
  filter(donor_id %in% clon_clus.sub$donor_id) %>%
  filter(donor_id %in% donor_md[donor_md$cohort!='Non',]$donor_id) %>%
  filter(!donor_id %in% donors.longit)

clustered = clustered %>%
  mutate(clustered = clone_id %in% clon_clus.sub$clone_id) %>%
  select(donor_id,clone_id,clustered) %>% distinct() %>%
  group_by(donor_id) %>% summarise(
    clones=n(), clusters=sum(clustered), clustered=sum(clustered)/n()*100) %>%
  left_join(clustered %>% group_by(donor_id) %>% tally(name = 'cells')) %>%
  left_join(donor_md %>% select(donor_id,cohort) %>% distinct())

tcr_analysis$clus[['pct']] = clustered

# Non-viral differential expression ############################################

tcr_analysis[['nonviral']] = list()

viral_hits = read.csv('../resources/viral_hits.csv')

library(zellkonverter)
library(scran)
sce = readH5AD('data/import/Tcell.h5ad')

# Filter assayed donors, have clone ID, are patient
sce = sce[,sce$donor_id %in% viral_hits$donor_id]
sce = sce[,(sce$clone_id!='')&(!is.na(sce$clone_id))]
sce=sce[,sce$donor_id %in% donor_md[donor_md$cohort!='Non',]$donor_id]
# Attach batch info to block on
colData(sce) = as_tibble(colData(sce)) %>%
  left_join(read.csv('../resources/Foster_2024_batch_id.csv')) %>% DataFrame()
# Covariate: viral annotation
sce$viral = sce$clone_id %in% viral_hits$clone_id

# Test for viral/non-viral
viral_deg_test = function(sce.test){
  res = scran::pairwiseTTests(
    sce.test,
    sce.test$viral,
    block=sce.test$batch_id,
    assay.type='X'
  )
  res = res$statistics[[1]] %>% data.frame() %>%
    mutate(gene=rownames(res$statistics[[1]])) %>% as_tibble() %>%
    arrange(-logFC) %>% mutate(annot = ifelse(logFC < 0,'viral','nonvir'))
  res
}

tcr_analysis[['nonviral']][['all']] = viral_deg_test(sce)
tcr_analysis[['nonviral']][['CD8.TemTerm']] = viral_deg_test(sce[,sce$pheno=='CD8.TemTerm'])

# Save output ##################################################################

save(tcr_analysis, file = 'data/tcr_analysis.RData')
