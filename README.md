# Code to reproduce _Tumour-intrinsic features shape T-cell differentiation through myeloma disease evolution_

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



