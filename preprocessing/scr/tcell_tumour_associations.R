# T cell-tumour association analysis ###########################################

library(tidyverse)
obs = read.csv('data/Tcell_obs.csv')
tcr = read.csv('data/import/Tcell-scTCR.csv')
load('data/comp.RData')
load('data/ords.RData')
load('data')
donor_md = read.csv('data/import/metadata-donor.csv')



#T cell skewing / premature skewing / abundance /
#clonality / effector signature / non-vir
#VS
#paraprotein/B2m (omitting Maura), aspirate CD138, % tumour stress



#nonvir vs tumour pathway corr
