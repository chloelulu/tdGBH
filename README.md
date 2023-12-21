
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tdGBH

Two-Dimensional Group Benjamini-Hochberg (tdGBH) Procedure is a multiple
testing procedure that offers enhanced power over traditional FDR
control methods when the data presents two-way grouping structure. It
reweights the p-values based on the informativeness of the respective
grouping directions. It requires a matrix of p-values as input with rows
and columns corresponding to the two grouping directions (e.g. genes by
cell types). For details of the proposed method, please refer to the
following paper:

Lu Yang, Pei Wang, Jun Chen. (2023) 2dGBH: Two-dimensional Group
Benjamini-Hochberg Procedure for False Discovery Rate Control in Two-Way
Multiple Testing of Genomic Data.

## Installation

You can install tdGBH as follows:

``` r
# install.packages("devtools")
devtools::install_github("chloelulu/tdGBH")
```

## Example

The input matrix contains the raw p-values from testing the associations
between the intake of 214 different nutrients and the abundance of 37
bacterial genera (Wu et al., 2021, Science). The output is the 2dGBH
FDR-adjusted p-values (or q-values).

``` r
library(tdGBH)
data(P)
p.adj <- tdGBH(P)
```
