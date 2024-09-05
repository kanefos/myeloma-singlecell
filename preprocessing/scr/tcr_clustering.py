# TCR clustering #######################################################################################################

import os
import numpy as np
import pandas as pd
pd.set_option('display.width', 600)
pd.set_option('display.max_columns', 20)
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from tcrdist.repertoire import TCRrep
from tcrdist.sample import _default_sampler
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.centers import calc_radii
from tcrdist.public import _neighbors_sparse_variable_radius, _neighbors_variable_radius
from tcrdist.public import TCRpublic
from tcrdist.ecdf import _plot_manuscript_ecdfs
from tcrdist.pgen import OlgaModel
olga_beta  = OlgaModel(chain_folder = "human_T_beta", recomb_type="VDJ")
olga_alpha = OlgaModel(chain_folder = "human_T_alpha", recomb_type="VJ")

# Import and prep for tcrdist3 #########################################################################################

tcr = pd.read_csv('data/import/Tcell-scTCR.csv')

# Remove missing
for i in ['CDR3aa','v_gene','j_gene']:
    tcr=tcr.loc[~tcr[i].isna(),]
    tcr = tcr.loc[tcr[i]!='',]

# Expanded (>1) only
tcr = tcr.loc[tcr.clone_size>1,]

# Prep for tcrdist
tcr['subject'] = tcr['donor_id']
tcr = tcr[['subject','chain','CDR3aa','v_gene','j_gene']]
for i in ['v_gene','j_gene']:
    tcr[i]=tcr[i].astype(str)+'*01'
    
# Alpha/beta data
tcr = { 'a':tcr.loc[tcr.chain=='TRA'],'b':tcr.loc[tcr.chain=='TRB'] }

# Prep columns for tcrdist3
for c in ['a','b']:
    column_map = {'CDR3aa':'cdr3_'+c+'_aa','v_gene':'v_'+c+'_gene','j_gene':'j_'+c+'_gene'}
    tcr[c].rename(columns=column_map, inplace=True)


# Run tcrdist3 over chains and donors ##################################################################################

clustering_output = {}

for c in ["a","b"]:
    for d in tcr[c].subject.unique():

        # Instance TCRrep
        tr = TCRrep(
            cell_df = tcr[c].loc[(tcr[c].chain=={'a':'TRA','b':'TRB'}[c])&(tcr[c].subject==d)].copy(),
            organism = 'human',
            chains = [{'a':'alpha','b':'beta'}[c]],
            compute_distances = False
        )
        # Distance
        tr.compute_distances()
        # TCRpublic class for reporting publicities
        #fixed radius 18
        tp = TCRpublic(tcrrep = tr, output_html_name = "tmp.html")
        # N to be public (in this instance, all)
        tp.query_str = 'nsubject > 0'
        #by calling, .report() an html report is made
        public = tp.report()

        # The non-redundant groups of quasipublic clones
        clus_df = public['clone_df']
        # Count CDR3 per cluster, filter for >1 (are actual clusters)
        clus_df['n_cdr3s'] = clus_df['cdr3s'].apply(len)
        clus_df = clus_df.loc[clus_df['n_cdr3s']>1,]
        # Prep for export
        clus_df = clus_df[['subject', 'chain', 'cdr3_'+c+'_aa', 'v_'+c+'_gene',
                           'j_'+c+'_gene', 'cdr3s', 'n_cdr3s']]
        column_map = {'cdr3_'+c+'_aa':'CDR3aa','v_'+c+'_gene':'v_gene','j_'+c+'_gene':'j_gene'}
        clus_df.rename(columns=column_map, inplace=True)

        clustering_output[d+'-'+c] = clus_df

        # Delete `tmp.html`
        os.remove("tmp.html")

pd.concat(clustering_output.values(), ignore_index=True
    ).to_csv('data/tcr_clustering.csv',index=None)
    
