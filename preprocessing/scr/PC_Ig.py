# Plasma cell Ig clonality calculation #################################################################################

import scanpy as sc
import pandas as pd
pd.set_option('display.width', 600)
pd.set_option('display.max_columns', 20)
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Import and pre-processing ############################################################################################

adata = sc.read('data/import/PC.h5ad')
md = pd.read_csv('data/import/metadata-donor.csv')

# Filter cancer only and samples without PC-excluding sort
adata = adata[adata.obs.donor_id.isin(md.loc[md.cohort.isin(['MGUS','SMM','MM'])].donor_id),:]
adata = adata[adata.obs.study_id!='Zheng_2021',]
adata = adata[~adata.obs.donor_id.isin(['Foster_2024.MM1','Foster_2024.MM2'])]

# Filter N PC > 50 cells
donor_PC = adata.obs.donor_id.value_counts()[adata.obs.donor_id.value_counts()>50].index
adata = adata[adata.obs.donor_id.isin( donor_PC ), :]

# Get Ig genes for each group
# Variable: light K/L variable; heavy variable
Ig_genes = {}
for g in ['IGKV',"IGLV",'IGHV']:
    Ig_genes[g] = adata.var_names[adata.var_names.str.startswith((g))].to_list()
# Heavy constant
Ig_genes['IGHC'] = ['IGHA2','IGHE','IGHG4','IGHG2','IGHGP','IGHA1','IGHG1','IGHG3','IGHD','IGHM']

# Subset adata to Ig genes only
Ig_genes_all = []
for g in Ig_genes.values():
    Ig_genes_all=Ig_genes_all+g
adata = adata[:,adata.var_names.isin(Ig_genes_all)]


# Per-cell, annotate name of most highly-expressed Ig gene #############################################################

def per_cell_annot_ig(adata_sub, donor_name):
    adata_sub = adata_sub.copy()

    # Create output df
    outputDF = pd.DataFrame(index=adata_sub.obs_names)
    outputDF['donor_id'] = donor_name

    # For each set of Ig genes
    for IgGene in ['IGKV', "IGLV", 'IGHV', 'IGHC']:
        # Subset to genes, get genes
        adata_sub_gene = adata_sub[:, [g for g in Ig_genes[IgGene] if g in adata_sub.var_names]]
        sub_genes = adata_sub_gene.var_names.to_list()
        # get dense matrix
        csr_mat = adata_sub_gene.X
        mat = csr_mat.todense()
        # note cells w 0 expression
        rowSums = mat.sum(axis=1)
        zeroExpr = rowSums == 0
        # note cells with non-zero expression in a single gene (thus cannot have secondary chains)
        exprOneGene = ((mat != 0).sum(1) == 1)
        # index of most highly expressed gene, then the gene at this index
        indexMax = np.argmax(mat, axis=1)
        exprMax = mat.max(axis=1)
        geneMax = [sub_genes[i] for i in [i[0] for i in np.array(indexMax)]]
        # for secondary chains, get index of second most highly expressed gene, then the gene
        indexSecond = np.array(
            np.argsort(mat.transpose(), axis=0)[-2]).transpose()  # had to be done in transposed form, idk
        geneSecond = [sub_genes[i] for i in indexSecond.reshape(1, len(indexSecond)).tolist()[0]]

        # get N counts and proportion of all Ig genes is from most expressed chain
        outputDF[IgGene + 'nMax'] = exprMax.A1
        outputDF[IgGene + 'propMax'] = (exprMax.A1 / rowSums.A1)

        # write to output DF
        outputDF[IgGene + '1'] = geneMax
        outputDF[IgGene + '2'] = geneSecond

        # filter for thise cells w known no 2ndary or no chains
        # give cells w zero expression '' value
        outputDF.loc[[i[0] for i in zeroExpr.tolist()], IgGene + str('1')] = ''
        outputDF.loc[[i[0] for i in zeroExpr.tolist()], IgGene + str('2')] = ''
        # give cells expressing a single gene a '' value in second gene
        outputDF.loc[[i[0] for i in exprOneGene.tolist()], IgGene + str('2')] = ''

    # Process as Ab data

    # Selecting light chain
    # remove if equal counts (useless/weird)
    outputDF.loc[outputDF['IGKVnMax'] != outputDF['IGLVnMax']]
    # set the primary/secondary lightV as the most abundant
    outputDF['lightV1'] = np.where(outputDF['IGKVnMax'] > outputDF['IGLVnMax'], outputDF['IGKV1'], outputDF['IGLV1'])
    outputDF['lightV2'] = np.where(outputDF['IGKVnMax'] > outputDF['IGLVnMax'], outputDF['IGKV2'], outputDF['IGLV2'])

    # and light max prop
    outputDF['lightV1prop'] = np.where(outputDF['IGKVnMax'] > outputDF['IGLVnMax'],
                                       outputDF['IGKVpropMax'], outputDF['IGLVpropMax'])

    return outputDF[['donor_id', 'lightV1prop', 'lightV1', 'lightV2', 'IGHVpropMax',
                     'IGHV1', 'IGHV2', 'IGHCpropMax', 'IGHC1', 'IGHC2']]

# Iterate over donors
donor_results = {}
for donor in adata.obs.donor_id.unique():
    donor_results[donor] = per_cell_annot_ig(adata[adata.obs.donor_id==donor,],donor)
out = pd.concat(donor_results.values(), ignore_index=False)

# Calculate the counts and proportion of each Ig gene per donor ########################################################

# Donor total cells
out_donor_ncells = out.donor_id.value_counts().reset_index().rename(columns={'index': 'donor_id', 'donor_id': 'donor_ncells'})

# Per-gene dict
donor_ig_gene_count_dict = {}

# For each donor, for each Ig gene, calc prop of gene in all cells, select the top (>33%)
# If the top has a mean per-cell prop <75% (suspect primary/second), select 2nd also
for ig_gene in ['lightV1', 'lightV2', 'IGHV1', 'IGHV2', 'IGHC1', 'IGHC2']:

    donor_ig_gene_count = out.groupby(['donor_id', ig_gene]).count(
        ).rename(columns={'lightV1prop': 'count'}).reset_index(
        ).loc[:, ['donor_id', ig_gene, 'count']].sort_values(by=["count"], ascending=False).merge(out_donor_ncells)

    # calc prop
    donor_ig_gene_count['prop'] = donor_ig_gene_count['count'] / donor_ig_gene_count['donor_ncells']
    # filter zeroes (unhelpful)
    donor_ig_gene_count = donor_ig_gene_count.loc[donor_ig_gene_count['count'] > 0, :]

    # set and rename ig gene
    donor_ig_gene_count['ig'] = ig_gene
    donor_ig_gene_count = donor_ig_gene_count.rename(columns={ig_gene: 'gene'})

    # output
    donor_ig_gene_count_dict[ig_gene] = donor_ig_gene_count

# Combine all
donor_ig_count = pd.concat(donor_ig_gene_count_dict.values(), ignore_index=False)

# Define a malignant Ig signature ######################################################################################

# Remove highly non-clonal (lowly-annotated) Ig
donor_ig_count_top = donor_ig_count.loc[donor_ig_count.prop > 0.1, :]

# Nested dict of donor/Ig_gene/annotated
donor_ig_top = {}
for donor_i in donor_ig_count_top.donor_id.unique():

    donor_ig_top[donor_i] = {}

    for ig_gene in donor_ig_count_top.ig.unique():

        donor_ig_gene_out = donor_ig_count_top.loc[
                            (donor_ig_count_top.donor_id == donor_i) & (donor_ig_count_top.ig == ig_gene), :][
                                'gene'].to_list()[:1]

        donor_ig_gene_out = list(g for g in donor_ig_gene_out if g != '')

        donor_ig_top[donor_i][ig_gene] = donor_ig_gene_out

# Identify putatative clonal plasma cells ##############################################################################

# Filter donors for whom Ig signatue derived
mPC_id = out.loc[out.donor_id.isin(donor_ig_top.keys()), :]

# Iterate over each cell, compare expressed Ig genes with summation of abundance in entire individual
for i, r in mPC_id.iterrows():
    donor = r['donor_id']

    ig_gene_bool = {}
    for ig in ['lightV', 'IGHV', 'IGHC']:

        ig1 = ig + str(1)
        ig2 = ig + str(2)

        # If primary or secondary chain are in either the top 1st/2nd most abundance genes for that individual
        # And this is true across multiple ig genes (if there is data)
        # Then, predict clonal

        if r[ig1] != '':
            ig1_bool = bool(r[ig1] in donor_ig_top[donor][ig1]) | bool(r[ig1] in donor_ig_top[donor][ig2])
        else:
            ig1_bool = False

        if r[ig1] != '':
            ig2_bool = bool(r[ig2] in donor_ig_top[donor][ig1]) | bool(r[ig2] in donor_ig_top[donor][ig2])
        else:
            ig2_bool = False

        mPC_id.loc[i, ig + '_bool'] = ig1_bool | ig2_bool

# Define clonal (putatative malignant) PC (mPC)
# Expression of clonal Ig in min one chain
mPC_id['mPC'] = (mPC_id['lightV_bool']==False) & (mPC_id['IGHV_bool']==False) & (mPC_id['IGHC_bool']==False)
mPC_id['mPC'] = ~mPC_id['mPC']

mPC_id.to_csv('data/mPC_id.csv')