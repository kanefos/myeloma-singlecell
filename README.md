1.Rmd# Code to reproduce _Tumour-intrinsic features shape T-cell differentiation through myeloma disease evolution_

As single-cell data from [Sklavenitis-Pistofidis et al.](https://doi.org/10.1016/j.ccell.2022.10.017) and extended clinical data from [Maura et al.](https://doi.org/10.1038/s43018-023-00657-1) (i.e. paraprotein) are omitted (and must be acquired through direct contact with these authors), this repository will not perfectly recreate figures seen in the manuscript. However, the underlying analysis and plotting code is unchanged, so acquiring and integrating this data will generate manuscript figures.

## Overview

### integration

_Overview of methods vs directly implementable, but could be_

Contents
- .yaml
- workbooks

Extra inputs: 
- `raw_data.h5ad`? Pre-processed/labelled identically to the Zenodo data, does not exist in dataset
- Validation data, raw and processed (processed is included in this dir)

Yields Zenodo .h5ad files

### preprocessing

_All being input onto Snakefile_

1. `import.R` downloads data from Zenodo to local environment
2. `process.py` generates working objects from .h5ad
3. `abundances.R` generates `comp` and `ords` objects

Other input data
- Granja CITE-seq
- External labels?
- TCR database (used version)
- CoMMpass _(or, omit?)_
- Another Zenodo for all this?


### notesbooks / figures

- Takes output of above, generates .Rmd and saved plots



