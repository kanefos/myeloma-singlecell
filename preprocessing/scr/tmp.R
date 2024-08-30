library(tidyverse)
valid = read.csv('../integration/data/validation_Tcell-obs.csv') %>% 
  left_join(read.csv('../integration/data/validation_Tcell.csv'))

head(valid,2)



# obs
obs = read.csv('data/validation_data/Tcell2_validT_run1-obs.csv') %>% 
  left_join(read.csv('data/validation_data/Tcell2_validT_run1.csv')) %>% 
  filter(study=='Botta2023') %>% 
  filter(!pheno_pred %in% c('Unknown','Myeloid','CD3_neg','Doublet')) 

# md 
md = as_tibble(read.csv('data/validation_data/validation-metadata.csv'))
valid_names_neat=data.frame(study=c('Botta2023','Lasry2023'),
                            name_neat=c('Botta et al. (2023)','Lasry et al (2023)'))



valid = read.csv('../integration/data/validation_obs.csv')


# comp
thr=0.2 #as recommended
valid_comp = read.csv(paste0('data/validation_data/comp-CertThr',thr,'.csv')) %>% as_tibble()
# ords
pca = read.csv(paste0('data/validation_data/pca-CertThr',thr,'.csv')) %>% as_tibble()

# scores
load(file='data/validation_data/validT_20240212-gsT.RData') #scr
