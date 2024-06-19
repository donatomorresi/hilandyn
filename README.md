
<!-- README.md is generated from README.Rmd. Please edit that file -->

# High-dimensional detection of Landscape Dynamics (HILANDYN)

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-1.0.5.9000-blue.svg)](https://github.com/donatomorresi/hilandyn)
[![R build
status](https://github.com/donatomorresi/hilandyn/workflows/R-CMD-check/badge.svg)](https://github.com/donatomorresi/hilandyn/actions)
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

Morresi, D., Maeng, H., Marzano, R., Lingua, E., Motta, R., & Garbarino,
M. (2024). High-dimensional detection of Landscape Dynamics: a Landsat
time series-based algorithm for forest disturbance mapping and beyond.
GIScience & Remote Sensing, 61(1), 2365001.

Maeng, H. (2019). Adaptive multiscale approaches to regression and trend
segmentation. Ph.D. thesis, London School of Economics and Political
Science.
