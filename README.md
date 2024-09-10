# Code to reproduce _Tumour-intrinsic features shape T-cell differentiation through myeloma disease evolution_

As single-cell data from [Sklavenitis-Pistofidis et al.](https://doi.org/10.1016/j.ccell.2022.10.017) and extended clinical data from [Maura et al.](https://doi.org/10.1038/s43018-023-00657-1) (i.e. paraprotein) are omitted (and must be acquired through direct contact with these authors), this repository will not perfectly recreate figures seen in the manuscript. However, the underlying analysis and plotting code is unchanged, so acquiring and integrating this data will generate manuscript figures.

## 1. `/integration`

- _Overview of methods vs directly implementable._
- `raw_data.h5ad` will be made available at publication or through communication with editor.

## 2. `/preprocessing`

Generates processed and analysed data for plotting from [Zenodo](https://doi.org/10.5281/zenodo.13646014).

Omitted input data
- Granja et al. CITE-seq data
- CoMMpass data (requires controlled access)

## 3. `/notesbooks`, `/figures`

One `/preprocessing` is run in entirety, each figure's notebook can be ran to view and save associated figures.
