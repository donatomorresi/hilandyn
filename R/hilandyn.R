#' hilandyn: Implementation of the High-dimensional detection of Landscape Dynamics (HILANDYN) algorithm for mapping forest disturbance dynamics through the temporal segmentation of high-dimensional Landsat time series. 
#' 
#' HILANDYN employs the High-dimensional Trend Segmentation (HiTS) procedure proposed by Maeng (2019) to detect changepoints in linear trends from Landsat time series including information in the spatial and spectral dimensions.
#' The HiTS procedure aims to detect changepoints in a piecewise linear signal where their number and location are unknown. Changes can occur in the intercept, slope or both of linear trends.
#' To start with, see the function \code{hilandyn_map}.
#' 
#' @author Donato Morresi \email{donato.morresi@@unito.it}, Hyeyoung Maeng \email{hyeyoung.maeng@@durham.ac.uk}, Matteo Garbarino \email{matteo.garbarino@@unito.it}
#' @references Maeng, H. (2019). Adaptive multiscale approaches to regression and trend segmentation. Ph.D. thesis, London School of Economics and Political Science.
#' @seealso \code{\link{hilandyn_map}}, \code{\link{hilandyn_int}}
#' @name hilandyn
#' @useDynLib hilandyn
NULL