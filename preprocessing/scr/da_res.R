# Differential abundance testing ###############################################

library(tidyverse)
source('../resources/aes.R')
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
load(file='data/comp.RData')
load(file='data/ords.RData')
results = list(panImm=list(),Tcell=list())

# Functions ####################################################################

# Cluster-level compositional linear mixed-effect model
celltype.LMM <- function(
    data.in,form.test,
    celltype='celltype',cond='cohort',randomEffect="donor_id",abundance='clr'
){
  # Iterate over input celltypes
  celltypes = data.in[,celltype] %>% pull() %>% unique()
  celltypes.res.list = list()
  for ( i in celltypes ) {

    data.in.sub = data.in[data.in[,celltype]==i,]

    # Skip if zero values
    if ( dim(data.in.sub)[1]==0 ){next}

    # The data defines the model fit

    # If N donor != N sample (pseudoreplication), fit linear mixed effect model
    if ( length(unique(data.in.sub$sample_id)) != length(unique(data.in.sub$donor_id)) ){
      # Re-define formula for randomEffect
      form.test.LMM = paste0(form.test," + (1|",randomEffect,")")
      # Try fitting LMM, if any errors (likely due to randomEffect) fit linear model
      model = tryCatch(
        {
          lmerTest::lmer(form.test.LMM, data=data.in.sub, REML=F)
        },
        error=function(e){
          lm(form.test, data=data.in.sub)
        }
      )
    # If not pesudoreplication (N donor = N sample), fit regular linear model
    } else if (length(unique(pull(data.in.sub[,cond])))>1) {
      model = lm(form.test, data=data.in.sub)
      # If anything (?) else just skip
    } else {
      next
    }

    # Calculate p-value via ANOVA
    anova.df = data.frame(anova(model))
    as_tibble(anova.df) %>%
      dplyr::rename('p'=rev(colnames(anova.df))[1]) %>%
      mutate(cond = rownames(anova.df)) %>% select(cond,p) %>%
      pivot_wider(names_from=cond, values_from=p) %>%
      setNames(paste0('p.', names(.))) -> pval

    # Calculate absolute difference in median CLR
    data.in.sub[,'cond'] = data.in.sub[,cond]
    data.in.sub %>%
      group_by(cond) %>%
      summarise(median=median(!!as.symbol(abundance))) %>% #, mean=mean(!!as.symbol(abundance))) %>%
      pivot_wider(names_from='cond',values_from='median') %>%
      mutate(diff=.[2]-.[1]) %>%
      pull(diff) -> diff

    # log2FC on 0-1 scaled CLR
    data.in.sub %>%
      group_by(cond) %>% mutate(clr = range01(clr)) %>%
      summarise(median=median(!!as.symbol(abundance))) %>% #, mean=mean(!!as.symbol(abundance))) %>%
      pivot_wider(names_from='cond',values_from='median') %>%
      mutate(x=log2(.[1]/.[2])) %>%
      .[3] %>% pull(.) -> log2fc

    # Create dataframe
    celltypes.res.list[[i]] = pval %>%
      mutate(diff=diff[1,1],log2fc=log2fc[1,1])
  }

  # Combine and output
  celltypes.res = bind_rows(celltypes.res.list, .id='celltype')
  return(celltypes.res)

}

# Cluster-level correlation
celltype.cor = function(data.in, cont.var=''){
  celltypes.in = unique(data.in$celltype)

  celltypes=c()
  pvals=c()
  cors=c()
  for ( celltype.i in celltypes.in ){
    data.in.celltype = data.in[data.in$celltype==celltype.i,]

    cor.res = stats::cor.test(pull(data.in.celltype,cont.var),data.in.celltype$clr)
    celltypes = c(celltypes,celltype.i)
    pvals = c(pvals,cor.res$p.value)
    cors = c(cors,cor.res$estimate)
  }
  as_tibble(data.frame(celltype=celltypes,pval=pvals,cor=cors))
}

# Across tissue  ###############################################################

results$panImm[['tissue']] = celltype.LMM(
  comp$panImm, 'clr ~ tissue',cond='tissue') %>% arrange(p.tissue)

results$Tcell[['tissue']] = celltype.LMM(
  comp$Tcell, 'clr ~ tissue',cond='tissue') %>% arrange(p.tissue)

# Across cohorts  ##############################################################

# Test several hypothesis
test.formulae = list(alone="clr ~ test.var",addAge="clr ~ test.var + age")
test.comparisons = list(
  Non.MGUS = c('Non','MGUS'), Non.SMM = c('Non','SMM'),
  Non.MM = c('Non','MM'), SMM.MM=c('SMM','MM'))

# Diagnosis, panImm
results$Tcell[['cohort']] = list()
dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100)
cohort.n = dat %>% select(donor_id,cohort) %>% distinct() %>% group_by(cohort) %>% tally()
for (formula in names(test.formulae)) {
  #print(formula)
  #print(test.formulae[[formula]])
  results$panImm[['cohort']][[formula]]=list()
  for (test.var in names(test.comparisons)){
    #print(test.var)
    #print(test.comparisons[[test.var]])
    dat.test = dat %>% filter(cohort %in% test.comparisons[[test.var]])
    dat.test[['cohort']] = factor(dat.test[['cohort']], levels= test.comparisons[[test.var]])

    res = celltype.LMM(dat.test, 'clr ~ cohort',cond='cohort') %>% arrange(p.cohort)
    res$denom = as.character(test.comparisons[[test.var]][1])
    res$denom.n = cohort.n[cohort.n$cohort==res$denom[1],]$n
    res$numer = as.character(test.comparisons[[test.var]][2])
    res$numer.n = cohort.n[cohort.n$cohort==res$numer[1],]$n

    results$panImm[['cohort']][[formula]][[test.var]] = res
  }
}

# Diagnosis, Tcell
results$Tcell[['cohort']] = list()
dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100)
cohort.n = dat %>% select(donor_id,cohort) %>% distinct() %>% group_by(cohort) %>% tally()
for (formula in names(test.formulae)) {
  #print(formula)
  #print(test.formulae[[formula]])
  results$Tcell[['cohort']][[formula]]=list()
  for (test.var in names(test.comparisons)){
    #print(test.var)
    #print(test.comparisons[[test.var]])
    dat.test = dat %>% filter(cohort %in% test.comparisons[[test.var]])
    dat.test[['cohort']] = factor(dat.test[['cohort']], levels= test.comparisons[[test.var]])

    res = celltype.LMM(dat.test, 'clr ~ cohort',cond='cohort') %>% arrange(p.cohort)
    res$denom = as.character(test.comparisons[[test.var]][1])
    res$denom.n = cohort.n[cohort.n$cohort==res$denom[1],]$n
    res$numer = as.character(test.comparisons[[test.var]][2])
    res$numer.n = cohort.n[cohort.n$cohort==res$numer[1],]$n

    results$Tcell[['cohort']][[formula]][[test.var]] = res
  }
}

# Diagnosis SMM/MM vs true healthy
test.comparisons = list(
  HD.SMM = c('HD','SMM'), HD.MM = c('HD','MM'))
dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
  filter(cohort %in% c('SMM','MM') | noncancer_sampling=='HD')
dat$cohort[dat$cohort=='Non']='HD'
cohort.n = dat %>% select(donor_id,cohort) %>% distinct() %>% group_by(cohort) %>% tally()
for (formula in names(test.formulae)) {
  #print(formula)
  #print(test.formulae[[formula]])
  for (test.var in names(test.comparisons)){
    #print(test.var)
    #print(test.comparisons[[test.var]])
    dat.test = dat %>% filter(cohort %in% test.comparisons[[test.var]])
    dat.test[['cohort']] = factor(dat.test[['cohort']], levels= test.comparisons[[test.var]])

    res = celltype.LMM(dat.test, 'clr ~ cohort',cond='cohort') %>% arrange(p.cohort)
    res$denom = as.character(test.comparisons[[test.var]][1])
    res$denom.n = cohort.n[cohort.n$cohort==res$denom[1],]$n
    res$numer = as.character(test.comparisons[[test.var]][2])
    res$numer.n = cohort.n[cohort.n$cohort==res$numer[1],]$n

    results$Tcell[['cohort']][[formula]][[test.var]] = res
  }
}

# With age #####################################################################

results$Tcell[['age']] = list()

dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
  filter(!is.na(age))

results$Tcell$age[['all']] = dat %>% celltype.cor(.,'age')
for (i in unique(dat$cohort)){
  results$Tcell$age[[i]] = dat %>% filter(cohort==i) %>% celltype.cor(.,'age')
}

# With paraprotein #############################################################

dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
  left_join(
    read_csv('data/import/metadata-donor.csv') %>%
      select(donor_id,Paraprotein)
  ) %>% filter(!is.na(Paraprotein))

results$Tcell[['Paraprotein']] = celltype.cor(dat,'Paraprotein') %>% arrange(pval)

# With aspirate % CD138+ cells #################################################

dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
  left_join(
    read_csv('data/import/metadata-donor.csv') %>%
      select(donor_id,aspirate_CD138pos)
  ) %>% filter(!is.na(aspirate_CD138pos), study_id!='Zavidij_2020')

results$Tcell[['aspirate_CD138pos']] = celltype.cor(dat,'aspirate_CD138pos') %>% arrange(pval)

# With tumour cell pathway expression ##########################################

scr = read_csv('data/tumour_modules_pct.csv')

dat = comp$Tcell %>% filter(
  !donor_id %in% c(donors.longit,donor.Tex_hi),tissue=='BM',sampleSize>100) %>%
  left_join(scr) %>% filter(!is.na(pct)) %>%
  select(sample_id,donor_id,celltype,clr,pathway_neat,pct)

scr_res = list()
for ( pw in unique(scr$pathway_neat)){
  scr_res[[pw]] = dat %>% filter(pathway_neat==pw) %>% celltype.cor(.,'pct')
}

results$Tcell[['tumour_modules_pct']] = bind_rows(scr_res,.id='pathway_neat')

# Save results #############################################################

save(results, file='data/da_results.RData')
