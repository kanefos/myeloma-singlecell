# Aesthetics
# ##########

#https://github.com/CartoDB/cartocolor

library(tidyverse)
theme_set(theme_classic())
library(ggpubr)
library(RColorBrewer)
library(viridis)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#Custom ggplot saving
#g = ggplot(mtcars, aes(mpg, drat, fill=as.character(gear), size=disp))+geom_point(shape=21)+th
#setwd('~/Desktop/')
#ggsave(plot=g, filename = "pl.svg",width=12, height=8, units='cm', bg = "transparent", dpi = 300)
#ggsave(plot=g, filename = "pl.pdf",width=12, height=8, units='cm', bg = "transparent", dpi = 300)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Selected donors
donors.longit = c('Foster_2024.MM4',"Foster_2024.MM3",'Foster_2024.SMM13',
                  'Liu_2021.37692_SMM','Liu_2021.47491_SMM',
                  'Liu_2021.58408_SMM','Oetjen_2018.C_T2','Oetjen_2018.S_T2')
donor.Tex_hi = 'Zheng_2021.P20190122'

# Cohort colors
disease_col = c('HD'='#77c8e0','HIP'='cadetblue4','MGUS'='#f2d37c','SMM'='#cf8a15','MM'='#eb676e')
diagnosis_col = c('Non'='#77c8e0','Non-cancer'='#77c8e0','HD'='#77c8e0','HIP'='cadetblue4','MGUS'='#f2d37c','SMM'='#cf8a15','MM'='#eb676e')
cohort_col = c('Non'='#77c8e0','Non-cancer'='#77c8e0','HD'='#77c8e0','HIP'='cadetblue4','MGUS'='#f2d37c','SMM'='#cf8a15','MM'='#eb676e')
smm_risk_col = c('high'='#eb7a09','low'='#e3a336')

# Study colors !!!!!
#study_colors = setNames(
#  c(colorRampPalette(brewer.pal(8, "Set2"))(12),colorRampPalette(brewer.pal(8, "Set1"))(2)),
#  read.csv('data/metadata/study_md.tsv',sep='\t')$study
#)
#study_colors['this_study']=study_colors['COSMOS21']
#study_name_neat_colors = setNames(
#  c(colorRampPalette(brewer.pal(8, "Set2"))(12),colorRampPalette(brewer.pal(8, "Set1"))(2)),
#  read.csv('data/metadata/study_md.tsv',sep='\t')$name_neat
#)

# Tissue colors
tissue_colors=c('BM'='#4C71A4','PB'='#E41A1C','BM & PB'='#47A265')

# pan-Immune celltype colors
bm_lineage_colors = setNames(
  c("#4C71A4","#E41A1C","#47A265","#7B7281","#CB6651","#FFAF13","#E8D430","#B05B3A","#F781BF"),
  c("T cell","Progenitor","PC","Myeloid","NK cell","B cell","Nonhaem","Platelet","Neutrophil")
)

# Tcell colors
cols_order = c("#7FC97F", "#FCCDE5", "#F0027F", "#BC80BD", "#FFAF13", "#FF7F00", "#FDC086",
               "#8DD3C7", "#A65628","#F781BF","#F781BF", "#377EB8", "#E41A1C", "#BEBADA", "#A6D854",
               "#CBD5E8", "#B3CDE3", "#4DAF4A", "#E5C494","#E5C494","#666666")
Tcell_ph_col = c("CD4.CTL","CD4.Tcm","Prolif.","CD8.TemActive","CD4.Tn","CD8.Tn","CD8.Tcm",
                "CD8.TEMRA","CD4.Th17","Treg","CD4.Treg","CD8.TemTerm","CD4.Tem","CD8.Tem.IL7R","CD8.Tem.KLRG1",
                "CD8.Trm","MAIT_gdT","ISG.ISG15","ISG.IFIT2",'Teff.IFIT2',"CD8.Tex")
Tcell_pheno_colors = set_names(cols_order,Tcell_ph_col)
Tcell_order = c(
  "CD4.Tn","CD4.Tcm","CD4.Th17","Treg",'CD4.Treg',"CD4.Tem","CD4.CTL",
  "CD8.Tn","CD8.Tcm","CD8.Tem.IL7R","CD8.Tem.KLRG1","CD8.TemActive","CD8.Trm","CD8.Tex","CD8.TEMRA","CD8.TemTerm",
  "MAIT_gdT","Prolif.","ISG.ISG15","ISG.IFIT2",'Teff.IFIT2'
)
Tcell_subsets = c('CD4+','CD4+','CD4+','CD4+','CD4+','CD4+','CD4+',
                 'CD8+','CD8+','CD8+','CD8+','CD8+','CD8+','CD8+','CD8+','CD8+',
                 'Other','Other','Other','Other','Other')
cell_order = rev(c(
  "CD4.Tn","CD4.Tn.ISG","CD4.Tcm","CD4.TcmActive","CD4.Th17","CD4.Tem","CD4.CTL","CD4.Treg",
  "CD8.Tn",'CD8.Tcm',"CD8.Tem.IL7R","CD8.Tem.CCR9","CD8.Tem.KLRG1","CD8.TemActive",
  "CD8.Tem.JUND","CD8.Trm","CD8.Tex","CD8.TemTerm","CD8.NKT.TEMRA","CD8.Tem.ISG",
  'MAIT','Prolif.','QC','Doublet'
))
pheno_subset_group = rev(c(
  "CD4+","CD4+","CD4+","CD4+","CD4+","CD4+","CD4+","CD4+",
  "CD8+",'CD8+',"CD8+","CD8+","CD8+","CD8+",
  "CD8+","CD8+","CD8+","CD8+","CD8+","CD8+",
  'Other','Other','Other','Other'
))

# mPC/hPC colors
mPC_col = c('True'='firebrick','False'='darkturquoise',`TRUE`='firebrick',`FALSE`='darkturquoise',
            'mPC'='firebrick','hPC'='darkturquoise', 'Clonal'='firebrick','Non-clonal'='darkturquoise')

# Long color vector
library(RColorBrewer)
qual_col_pals = RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual',]
long_col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Study colors
study_colors = setNames(
  rev(long_col_vector)[1:11],
  c("Foster_2024","Oetjen_2018","Bailur_2019","Zavidij_2020","Kfoury_2021",
           "Granja_2019","Zheng_2021","Liu_2021","Maura_2023","Conde_2022","Stephenson_2021"))

# Sort colors
sort_colors = setNames(
  long_col_vector[1:7],
  c("T cell-enriched/depleted","Unsorted","CD138-","CD138-CD45+",
    "CD235-","T cell-enriched","CD8-enriched")
)

# Custom ggplot
library("scales")
reverselog_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}
