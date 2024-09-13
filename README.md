## Code to reproduce _Tumour-intrinsic features shape T-cell differentiation through myeloma disease evolution_

Single-cell data from [Sklavenitis-Pistofidis et al.](https://doi.org/10.1016/j.ccell.2022.10.017) and extended clinical data from [Maura et al.](https://doi.org/10.1038/s43018-023-00657-1) (i.e. paraprotein) are omitted (and must be acquired through direct contact with these authors), this repository will not perfectly recreate figures seen in the manuscript. However, the underlying analysis and plotting code is unchanged, so acquiring and integrating this data will generate manuscript figures.

### 1. `/environments`

- `Rpackages.html`, install correct R package versions.
- `.yml` files used to construct conda env for analysis:
```shell
conda env create --name integration_env --file=integration_env.yml
conda env create --name preprocessing --file=preprocessing.yml
conda env create --name tcrdist3_env --file=tcrdist3_env.yml
```

### 2. `/integration`

- _Overview of methods vs directly implementable._
- `/environments/integration_env.yml` can be used to construct conda environment used for integration.
- `raw_data.h5ad` will be made available at publication or (for reviewers) through communication with editor.

### 3. `/preprocessing`

Generates processed and analysed data for constructing figures using data available on [Zenodo](https://doi.org/10.5281/zenodo.13646014).

Run:
```shell
# First, edit `run_preprocessing.py` to user python (from above conda env) and R locations
# Then, in bash shell, run within /myeloma-singlecell/preprocessing
python3 run_preprocessing.py
```

Omitted data:
- Single-cell data from [Sklavenitis-Pistofidis et al.](https://doi.org/10.1016/j.ccell.2022.10.017).
- Extended clinical data from [Maura et al.](https://doi.org/10.1038/s43018-023-00657-1).
- Granja et al. CITE-seq data, [available here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369).
- CoMMpass data (requires controlled access)

### 4. `/notebooks` & `/figures`

Generate all figures
```shell
# bash shell, within /myeloma-singlecell
Rscript make_figures.R
```





