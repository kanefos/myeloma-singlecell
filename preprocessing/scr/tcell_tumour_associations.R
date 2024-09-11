# T cell-tumour association analysis ###########################################

library(tidyverse)
source('../resources/aes.R')
obs = read.csv('data/Tcell_obs.csv')
tcr = read.csv('data/import/Tcell-scTCR.csv')
load('data/comp.RData')
load('data/ords.RData')
load('data/tcr_analysis.RData')
#load('data/tumour_modules_scr.RData')
tumour_modules_pct = read_csv('data/tumour_modules_pct.csv')
donor_md = read_csv('data/import/metadata-donor.csv')
sort_id = read.csv('../resources/Foster_2024_sort_id.csv')

# T cell-tumour feature-by-feature association #################################

# Data prep

# Mean Non-viral score
Nonviral_score = as_tibble(obs) %>% select(sample_id,donor_id,Nonviral) %>% drop_na() %>%
  filter(sample_id %in% (
    comp$Tcell %>% filter(sampleSize>100, !donor_id %in% donors.remove, tissue=='BM') %>%
      filter( ! sample_id %in% sort_id[sort_id$sort_id!='T cell-enriched',]$sample_id ) %>% pull(sample_id)
  )) %>% group_by(donor_id) %>% summarise(Nonviral=median(Nonviral))

# T effector signature
Teff_score = read_csv('data/Tcell_exp.csv') %>%
  left_join(obs %>% select(obs_names,donor_id)) %>%
  select(donor_id,CD69,IFNG,TNF) %>%
  group_by(donor_id) %>% summarise_all(.funs = mean) %>%
  mutate_at(c('CD69','IFNG','TNF'), funs(c(range01(.)))) %>%
  rowwise() %>% mutate(Teff_score = TNF+IFNG+CD69 / 3) %>%
  select(donor_id,Teff_score)

# Combining
donor_md %>%
  left_join(comp$Tcell %>% select(donor_id,study_id) %>% distinct()) %>%
  # Clinical data
  select(study_id,donor_id,Paraprotein,B2m,aspirate_CD138pos,cohort) %>%
  mutate(aspirate_CD138pos=ifelse(study_id=='Zavidij_2020',NA,aspirate_CD138pos)) %>%
  # Cancer only, non-longit
  filter(cohort!='Non', !donor_id %in% donors.longit) %>%
  # Tumour stress
  left_join(
    tumour_modules_pct %>% filter(p=='MP5.Stress') %>%
      mutate(pct.stress=pct) %>% select(donor_id,pct.stress)
  ) %>%
  # Clonality
  left_join(
    tcr_analysis$clonality %>% filter(group=='all') %>% select(donor_id,clon) %>% distinct()
  ) %>%
  # Teff score
  left_join((Teff_score)) %>%
  # T cell skewing and exagerated skewing
  left_join(
    ords$exag_skewing %>% select(sample_id,PC1,residuals) %>%
      left_join(comp$Tcell %>% select(sample_id,donor_id) %>% distinct()) %>%
      mutate(Tcell_skewing=PC1,exag_Tcell_skewing=residuals) %>%
      select(donor_id,Tcell_skewing,exag_Tcell_skewing)
  ) %>%
  # Non-viral score
  left_join(Nonviral_score) %>%
  # Cluster abundance
  left_join(
    comp$Tcell %>%
      filter(cohort!='Non',sampleSize>100, !donor_id %in% donors.remove, tissue=='BM') %>%
      filter( ! sample_id %in% sort_id[sort_id$sort_id!='T cell-enriched',]$sample_id ) %>%
      filter(celltype %in% c('ISG.IFIT2','CD8.TEMRA','CD8.TemTerm','CD8.Tn')) %>%
      select(donor_id,celltype,clr) %>%
      pivot_wider(names_from=celltype,values_from=clr)
  ) -> tcell_tum_merged

# Correlation testing
test_over = tcell_tum_merged %>% select(-study_id,-donor_id,-cohort) %>% colnames()
out = list()
for (i in test_over){
  for (j in test_over){
    if (i==j){next}
    x=tcell_tum_merged[[i]]
    y=tcell_tum_merged[[j]]
    xy = data.frame(x=x,y=y) %>% drop_na()
    out[[paste0(i,'_',j)]] = data.frame(
      i=i, j=j, p=cor.test(xy$x, xy$y)$p.value,r=cor.test(xy$x, xy$y)$estimate, n=dim(xy)[1])
  }
}

# Plotting prep
bind_rows(out) %>% as_tibble()
ft.tum = c('Paraprotein','B2m','aspirate_CD138pos','pct.stress')
ft.t = c('Tcell_skewing', 'exag_Tcell_skewing','CD8.TEMRA', 'clon',
         'Teff_score', 'Teff.IFIT2', 'Nonviral')
as_tibble(bind_rows(out)) %>% distinct() %>%
  mutate(fdr=p.adjust(p,method='BH')) %>%
  filter(i %in% ft.tum, j %in% ft.t) %>%
  # Rename neat
  mutate(i.pl=i,j.pl=j) %>%
  mutate(i.pl=replace(i.pl,i=='aspirate_CD138pos','Aspirate % CD138+')) %>%
  mutate(i.pl=replace(i.pl,i=='pct.stress','% tumour expressing\nstress pathway')) %>%
  mutate(j.pl=replace(j.pl,j=='clon',"Repertoire clonality" )) %>%
  mutate(j.pl=replace(j.pl,j=='CD8.TEMRA',"CD8.TEMRA\n(normalised % T cells)" )) %>%
  mutate(j.pl=replace(j.pl,j=='Teff.IFIT2',"Teff.IFIT2\n(normalised % T cells)" )) %>%
  mutate(j.pl=replace(j.pl,j=='Tcell_skewing',"T cell skewing" )) %>%
  mutate(j.pl=replace(j.pl,j=='exag_Tcell_skewing',"Exaggerated T cell aging" )) %>%
  mutate(j.pl=replace(j.pl,j=='Nonviral',"Non-viral specificity\nscore" )) %>%
  mutate(j.pl=replace(j.pl,j=='Teff_score',"Effector signature\n(CD69, TNF, IFNG)" )) %>%
  # Significance marks
  mutate(sig_gr = case_when(
    p >= 0.1 ~ "",
    p >= 0.05 & p < 0.1 ~ "",
    p >= 0.01 & p < 0.05 ~ "*",
    p >= 0.001 & p < 0.01 ~ "**",
    p < 0.001 ~ "***"
  )) -> tcell_tum_res

tcell_tum_res %>% write.csv('data/tcell_tum_association_results.csv',row.names=F)


# PC1 T cell skewing and tumour pathway expression #############################

dat = tumour_modules_pct %>% left_join(
    ords$Tcell$PCA %>% left_join(comp$Tcell %>% select(sample_id,donor_id) %>% distinct())
  ) %>%
  filter(!is.na(PC1))

out=list()
for (pw in unique(dat$pathway_neat)){
  dat2 = dat %>% filter(pathway_neat==pw)
  res = cor.test(dat2$pct, dat2$PC1)
  out[[pw]] = data.frame(p=res$p.value, r=res$estimate)
}
as_tibble(bind_rows(out, .id='pathway')) %>% mutate(fdr=p.adjust(p,'BH')) %>%
  write.csv('data/PC1_tumour_modules.csv',row.names=F)

# Non-viral T cell score and tumour pathway expression #########################

# Mean Non-viral score
Nonviral_score = as_tibble(obs) %>% select(sample_id,donor_id,Nonviral) %>% drop_na() %>%
  filter(sample_id %in% (
    comp$Tcell %>% filter(sampleSize>100, !donor_id %in% donors.remove, tissue=='BM') %>%
      filter( ! sample_id %in% sort_id[sort_id$sort_id!='T cell-enriched',]$sample_id ) %>% pull(sample_id)
  )) %>% group_by(donor_id) %>% summarise(Nonviral=median(Nonviral))

dat = tumour_modules_pct %>%
  left_join(Nonviral_score) %>% filter(!is.na(Nonviral))

out=list()
for (pw in unique(dat$pathway_neat)){
  dat2 = dat %>% filter(pathway_neat==pw)
  res = cor.test(dat2$pct, dat2$Nonviral)
  out[[pw]] = data.frame(p=res$p.value, r=res$estimate)
}
as_tibble(bind_rows(out, .id='pathway')) %>% mutate(fdr=p.adjust(p,'BH')) %>%
  write.csv('data/Nonviral_tumour_modules.csv',row.names=F)


# Paraprotein and T cell cluster abundance #####################################

IgG_and_IgA_isotype_tumours = read.csv('../resources/IgG_and_IgA_isotype_tumours.csv')$donor_id
Paraprotein_plot = donor_md %>% select(donor_id,Paraprotein) %>% drop_na() %>%
  filter(donor_id %in% IgG_and_IgA_isotype_tumours)

out=list()
for (i in unique(comp$Tcell$celltype)){

  dat = comp$Tcell %>% filter(
    !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
    filter( ! sample_id %in% sort_id[sort_id$sort_id!='T cell-enriched',]$sample_id ) %>%
    filter(celltype==i) %>% select(-Paraprotein) %>%
    left_join(Paraprotein_plot) %>% filter(!is.na(Paraprotein))

  res = cor.test(dat$clr, dat$Paraprotein)
  out[[i]] = data.frame(p=res$p.value, r=res$estimate)
}
as_tibble(bind_rows(out, .id='celltype')) %>% mutate(fdr=p.adjust(p,'BH')) %>%
  write.csv('data/Paraprotein_Tcell_correlation.csv',row.names=F)

# Tumour pathway expression and T cell cluster abundance #######################

out=list()
for (pw in unique(tumour_modules_pct$pathway_neat)){
out[[pw]]=list()
for (i in unique(comp$Tcell$celltype)){

dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
  filter( ! sample_id %in% sort_id[sort_id$sort_id!='T cell-enriched',]$sample_id ) %>%
  filter(celltype==i) %>%
  left_join(
    tumour_modules_pct %>% filter(pathway_neat==pw) %>% select(donor_id,pct)
  ) %>% filter(!is.na(pct))

  res = cor.test(dat$clr, dat$pct)
  out[[pw]][[i]] = data.frame(p=res$p.value, r=res$estimate)
}
out[[pw]]=bind_rows(out[[pw]],.id='celltype')
}

as_tibble(bind_rows(out, .id='pathway_neat')) %>% mutate(fdr=p.adjust(p,'BH')) %>%
  write.csv('data/Tcell_abundance_tumour_modules.csv',row.names=F)










