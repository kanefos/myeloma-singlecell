import os
import pandas as pd
import scanpy as sc

# Process .h5ad data ###########################################################

markers = pd.read_csv('../resources/markers.csv')

#panImm_obs
#PC_subset
#panImm_exp 
adata = sc.read('data/import/panImmune.h5ad')
adata.obs.to_csv('data/panImm_obs.csv',index=False)
adata[adata.obs.lineage=='PC',].write('data/import/PC.h5ad',compression='gzip')
sc.get.obs_df(adata, markers.loc[markers['plot']=='panImm',]['gene'].to_list()
              ).to_csv('data/panImm_exp.csv',index=False)

#Tcell_obs
#Tcell_exp
adata = sc.read('data/import/Tcell.h5ad')
adata.obs.to_csv('data/Tcell_obs.csv',index=False)
sc.get.obs_df(adata, markers.loc[markers['plot']=='Tcell',]['gene'].to_list()
              ).to_csv('data/Tcell_exp.csv',index=False)
