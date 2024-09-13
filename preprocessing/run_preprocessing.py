import os

conda = '/Users/username/miniconda3/envs/'
R_loc = '' #keep blank if using default R environment

preprocessing='preprocessing'
tcrdist_env="tcrdist3_env"

print('import.R #########################################################')
os.system(R_loc+"Rscript scr/import.R")

print('process.py #########################################################')
os.system(conda+preprocessing+"/bin/python scr/process.py")

print('abundances.R #########################################################')
os.system(R_loc+"Rscript scr/abundances.R")

print('PC_Ig.py #########################################################')
os.system(conda+preprocessing+"/bin/python scr/PC_Ig.py")

print('PC_analysis.R #########################################################')
os.system(R_loc+"Rscript scr/PC_analysis.R")

print('da_res.R #########################################################')
os.system(R_loc+"Rscript scr/da_res.R")

print('tcr_clustering.py #########################################################')
os.system(conda+tcrdist_env+"/bin/python3.8 scr/tcr_clustering.py")

print('tcr_analysis.R #########################################################')
os.system(R_loc+"Rscript scr/tcr_analysis.R")

print('tcell_tumour_associations.R #########################################################')
os.system(R_loc+"Rscript scr/tcell_tumour_associations.R")
