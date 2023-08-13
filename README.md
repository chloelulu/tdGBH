
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tdGBH

The goal of tdGBH is to perform two-dimensional group Benjamini-Hochberg
procedure for false discovery rate control in Two-Way multiple testing.

## Installation

You can install tdGBH like so:

``` r
devtools::install_github("chloelulu/tdGBH")
```

## Example

The input matrix includes raw p-values that represent the association
between the intake of 214 different nutrients and the abundance of 37
bacterial genera. Results returns the tdGBH adjusted p values storaged
in a matrix with the same dimension as the input matrix.

``` r
# library(tdGBH)
# data(P)
# res <- tdGBH(P)
```

## Reference

Lu Yang, Jun Chen. (2023) 2dGBH: Two-dimensional Group
Benjamini-Hochberg Procedure for False Discovery Rate Control in Two-Way
Multiple Testing.
