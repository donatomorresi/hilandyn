
<!-- README.md is generated from README.Rmd. Please edit that file -->

# High-dimensional detection of Landscape Dynamics (HILANDYN) <img align="right" width="250" src="man/figures/logo.png">

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-1.0.6.9000-blue.svg)](https://github.com/donatomorresi/hilandyn)
[![R build
status](https://github.com/donatomorresi/hilandyn/workflows/R-CMD-check/badge.svg)](https://github.com/donatomorresi/hilandyn/actions)
<!-- badges: end -->

`hilandyn` is an *R* package that provides an implementation of the
High-dimensional detection of Landscape Dynamics algorithm (Morresi *et
al.* 2024) for mapping forest disturbance dynamics by segmenting
high-dimensional Landsat time series into linear trends.

High-dimensional Landsat time series include information from the
spatial and spectral dimensions and are analysed using a modified
version of the High-dimensional Trend Segmentation (HiTS) procedure
proposed by Maeng (2019). The HiTS procedure aims to detect changepoints
in a piecewise linear signal where their number and location are
unknown. Changes can occur in the intercept, slope or both of linear
trends.

`hilandyn` uses the `terra` package for raster and vector data
management and the `future` package for parallel computation.

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
