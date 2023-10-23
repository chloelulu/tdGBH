
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tdGBH

tdGBH is a multiple testing procedure that offers enhanced power over
traditional FDR control methods like BH or ST. This is achieved by
assigning data-informed weights to each hypothesis. Specifically, tdGBH
conducts a two-dimensional group Benjamini-Hochberg procedure, targeting
false discovery rate control within Two-Way multiple testing scenarios.
It requires a matrix of p-values as input, structured with features as
rows and outcomes as columns. For a comprehensive understanding of
tdGBH, please refer to the subsequent paper:

Lu Yang, Pei Wang, Jun Chen. (2023) 2dGBH: Two-dimensional Group
Benjamini-Hochberg Procedure for False Discovery Rate Control in Two-Way
Multiple Testing.

## Installation

You can install tdGBH as follows:

``` r
devtools::install_github("chloelulu/tdGBH")
```

## Example

The input matrix includes raw p-values that represent the association
between the intake of 214 different nutrients and the abundance of 37
bacterial genera. The results provide the 2dGBH adjusted p-values in a
matrix that mirrors the dimensions of the original input matrix.

``` r
library(tdGBH)
data(P)
p.adj <- tdGBH(P)
```
