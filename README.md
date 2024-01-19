
<!-- README.md is generated from README.Rmd. Please edit that file -->

# High-dimensional detection of Landscape Dynamics (HILANDYN)

<!-- badges: start -->
<!-- badges: end -->

## About

The `hilandyn` *R* package provides an implementation of the
High-dimensional detection of Landscape Dynamics (HILANDYN) algorithm
for mapping forest disturbance dynamics through the segmentation of
high-dimensional Landsat time series. High-dimensional Landsat time
series include information from the spatial and spectral dimensions and
are analysed using the High-dimensional Trend Segmentation (HiTS)
procedure proposed by Maeng (2019). The HiTS procedure aims to detect
changepoints in a piecewise linear signal where their number and
location are unknown. Changes can occur in the intercept, slope or both
of linear trends. `hilandyn` uses the `terra` package for raster and
vector data management and the `future` package for parallel
computation.

## Installation

You can install the development version of hilandyn from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("donatomorresi/hilandyn")
```

### References

Maeng, H. (2019). Adaptive multiscale approaches to regression and trend
segmentation. Ph.D. thesis, London School of Economics and Political
Science.
