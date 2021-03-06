# CCPLS <img src="https://user-images.githubusercontent.com/49115350/114253008-345bf600-99e3-11eb-8b62-4af53ed665f9.png" width="30%" align="right" />

R package for estimating cell-cell communications from spatial transcriptome data with single-cell resolution.

Please see [preprint of CCPLS](https://www.biorxiv.org/content/10.1101/2022.01.12.476034v1) for details.

## Getting started

### Install

Requirements: R version 4.2.0 or higher.

```
# Install dependent packages
install.packages(c("cluster", "circlize", "dplyr", "pls", "purrr", "Seurat", "stringr"), dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
# Install CCPLS
install.packages("devtools", dependencies = TRUE)
devtools::install_github("bioinfo-tsukuba/CCPLS")
 ```

### Quick exmaple

#### 1. Prepare arguments for CCPLS

```
# Here is demonstraion by included dataset in CCPLS package.
# Please prepare dataset and output directory for your purpose.
load(system.file("extdata", "dataset.Rdata", package = "CCPLS"))
output_dir <- "~/CCPLS_test"

# Note that 400 cells were randomly extracted by Seq-Scope data (Cho et al., 2021), and this demonstration cannot be interpreted biologically.
```

#### 2. Run CCPLS

```
result_CCPLS <- CCPLS::cellCellReg(exp_mat, coord_mat, annot_mat, output_dir)
```

#### 3. View reports by CCPLS

<img src="https://user-images.githubusercontent.com/49115350/148733504-73c78ba4-b8d1-4c31-9026-925c827ae5cb.png" width="50%"><img src="https://user-images.githubusercontent.com/49115350/148733531-6943ef3e-ba43-466b-8177-c8cfb7f2e166.png" width="50%">

## License
Copyright (c) 2022 Takaho Tsuchiya and Bioinformatics Laboratory, Faculty of Medicne, University of Tsukuba released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).
