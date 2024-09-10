# ENVS
#conda snakemake = for running snakemake
#Rscript ... = default R env, which runs normally. need to package up
#     collate all R packages/versions?...   


#RUN ENTIRE PREPROCESSING DIR
#
#
#


rule notesbooks_and_figures:
    input:
    "FIG-x.Rmd"
#repeat for all
output:
    "FIG-x.html"
directory("figures/FIG-x")
#repeat for all figs
R:
    "rmarkdown::render('notebooks/tmp.Rmd',params=list(args = 'myarg'))"
#this requires some fidelling to run on cmd line, but is possible:
"Rscript -e 'rmarkdown:: ...etc'"


rule consolidate_data:
    input:
    "ALL PREVIOUS GENERATED FILES? OR JUST THOSE PROCESSED?"
output:
    "dummy file to include in rule all"
"i.e. check all others installed?"
shell:
    "dummy line"
