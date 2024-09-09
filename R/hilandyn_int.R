#' High-dimensional detection of Landscape Dynamics engine
#'
#' This function is intended to analyse the data extracted from the cells within a single spatial kernel.
#'
#' Typically, \code{hilandyn_int()} is not directly called by the user.
#' It processes inter-annual high-dimensional Landsat time series to detect changes in spectral trends within a spatial kernel.
#' Changes in the intercept, slope or both of linear trends are detected using the High-dimensional Trend Segmentation (HiTS) procedure proposed by \insertCite{maeng2019adaptive;textual}{hilandyn}.
#' The HiTS procedure aims at detecting changepoints in a piecewise linear signal where their number and location are unknown.
#' High-dimensional time series can include single or multiple reflectance bands and spectral indices, hereafter referred to as bands.
#' Input data can be a \code{numeric} \code{vector} extracted from a \code{SpatRaster} or a \code{matrix}. Each row in the matrix corresponds to a variable of the high-dimensional time series.
#' Impulsive noise, i.e. outliers in the time series, can be removed through an iterative procedure by setting the \code{noise_iter_max} and \code{nob_init_min} parameters.
#' One-year gaps in the time series are filled using either linear interpolation or extrapolation.
#'
#' @param x numeric vector or matrix. Input time series where each row of the matrix contains time-ordered data relative to a band and a cell in the spatial kernel. Vectors are transformed into matrices.
#' @param nob numeric vector or matrix. Number of clear observations available per-cell (in rows) and per-year (in columns) for producing reflectance composites. Vectors are transformed into matrices.
#' @param nb integer. Number of input bands in the time series.
#' @param nc integer. The number of cells within the spatial kernel.
#' @param ny integer. Number of time steps (years) in the time series.
#' @param nr integer. Number of rows (variables) in the high-dimensional time series.
#' @param yrs numeric vector. The sequence of years to be analysed.
#' @param foc_ind numeric vector. Row indices of the bands in the focal cell.
#' @param cell_weights logical. Whether to compute cell-based weights within the spatial kernel.
#' @param ev logical vector. A vector containing \code{NA} values to be used in case of processing failure.
#' @param cng_dir numeric vector. Direction (either 1 or -1) of the spectral change caused by a disturbance in each band. It is computed as the difference between pre-disturbance and post-disturbance values.
#' @param th_const numeric. Constant value controlling the sensitivity to deviations from linearity in the HiTS procedure. Typical values are comprised in the interval \eqn{[0.7, 1.3]}.
#' @param noise_iter_max integer. The maximum number of iterations allowed for removing impulsive noise with the noise filter. The noise filter is disabled when the value is 0.
#' @param nob_init_min integer. The minimum number of clear observations available per time step in the first two time steps of the time series.
#' @param rmse logical. Determines whether to compute the root mean square error for each band in the focal cell. Original and estimated values are processed before noise filtering and gap filling.
#' @param use_last logical. Determines whether changepoints detected at the last time point are ignored or not.
#' @param as_list logical. If \code{TRUE} the output is a \code{list}. Otherwise, the output is a \code{numeric} \code{vector}.
#'
#' @return If \code{as_list} is \code{TRUE}, a \code{list} containing the following elements. Otherwise the result is a \code{numeric} \code{vector}.
#'   \item{D_MAX_YR}{Year corresponding to the disturbance event with the maximum spectral change magnitude \code{D_MAX_MG} (single value).}
#'   \item{D_MAX_MG}{Maximum spectral change magnitude of any disturbance computed as the median of bands and cells within the spatial kernel (single value).}
#'   \item{D_FRS_YR}{Year corresponding to the first disturbance event (single value).}
#'   \item{D_FRS_MG}{Spectral change magnitude of the first disturbance event computed as the median of bands and cells within the spatial kernel (single value).}
#'   \item{D_LST_YR}{Year corresponding to the last disturbance event (single value).}
#'   \item{D_LST_MG}{Spectral change magnitude of the last disturbance event computed as the median of bands and cells within the spatial kernel (single value).}
#'   \item{D_NUM}{Total number of disturbance events (single value).}
#'   \item{G_MAX_YR}{Year corresponding to the greening event with the maximum spectral change magnitude \code{G_MAX_MG} (single value).}
#'   \item{G_MAX_MG}{Maximum spectral change magnitude of greening computed as the median of bands and cells within the spatial kernel (single value).}
#'   \item{G_FRS_YR}{Year corresponding to the first greening event (single value).}
#'   \item{G_FRS_MG}{Spectral change magnitude of the first greening event computed as the median of bands and cells within the spatial kernel (single value).}
#'   \item{G_LST_YR}{Year corresponding to the last disturbance event (single value).}
#'   \item{G_LST_MG}{Spectral change magnitude of the last disturbance event computed as the median of bands and cells within the spatial kernel (single value).}
#'   \item{G_NUM}{Total number of greening events (single value).}
#'   \item{N_GAP}{Number of gaps in the time series, if any (single value).}
#'   \item{N_NOISE}{Number of years containing impulsive noise, if any (single value).}
#'   \item{RMSE}{Root mean square error of each band in the focal cell (one value per band).}
#'   \item{LEN}{Length of segments (one value per year).}
#'   \item{CPT_ID}{Type of change (one value per year). One of the following values: 101 (disturbance); 102 (greening); 9 (other change). Otherwise 0 for no change.}
#'   \item{CPT_FOC}{Number of changepoints detected in the focal cell. The minimum value is zero and the maximum corresponds to the number of bands.}
#'   \item{NOISE}{Position of the impulsive noise in the time series (one value per year).}
#'   \item{MED_REL_MAG}{Median relative magnitude computed using all the cells in the spatial kernel and the bands.}
#'   \item{EST}{Estimated values of the bands in the focal cell (one value per year and band).}
#'   \item{SLO}{Slope of the linear segments of the bands in the focal cell (one value per year and band).}
#'   \item{MAG}{Magnitude in absolute terms of the bands in the focal cell (one value per year and band).}
#'   \item{REL_MAG}{Magnitude in relative terms of the bands in the focal cell (one value per year and band).}
#'
#' @author Donato Morresi, \email{donato.morresi@@unito.it}
#'
#' @references
#' \insertRef{maeng2019adaptive}{hilandyn}
#'
#' @examples
#' library(hilandyn)
#'
#' # Load raster data
#' data(lnd_si)
#' lnd_si <- terra::rast(lnd_si)
#'
#' # Extract values from spatial kernel
#' cells <- terra::adjacent(lnd_si, 2000, directions = "8", include = TRUE)
#' v <- c(as.matrix(lnd_si[sort(cells)]))
#' 
#' # Set parametrs
#' bands <- c("MSI", "TCW", "TCA")
#' nc <- 9
#' years <- 1985:2020
#' foc_ind <- c(5, 14, 23)
#' cng_dir <- c(-1, 1, 1)
#'
#' # Process data
#' out <- hilandyn_int(v,
#'                     nob = rep(999, nc * length(years)),
#'                     nb = length(bands),
#'                     nc = nc,
#'                     ny = length(years),
#'                     nr = length(bands) * nc,
#'                     yrs = years,
#'                     foc_ind = foc_ind,
#'                     cng_dir = rep(cng_dir, each = nc),
#'                     as_list = TRUE)
#'                 
#' # Plot results
#' dim(v) <- c(27, length(years))
#' par(mfrow=c(length(bands), 1))
#' 
#' for(i in seq_along(bands)) {
#'   plot(v[foc_ind[i],], type = "l", col = 2, lwd = 2, ylab = "", xlab = "", xaxt = "n")
#'   lines(out$EST[i,], type = "l", col = 4, lwd = 2, lty = 1)
#'   abline(v = which(out$CPT_ID %in% c(101)), lty = 2, col = 2, lwd = 2)
#'   axis(1, at = seq_along(years), labels = years)
#'   title(main = paste(bands[i], "(RMSE =", out$RMSE[i], ")"), adj = 0)
#' }
#'
#' @export

hilandyn_int <- function(x, nob, nb, nc, ny, nr, yrs, foc_ind, cell_weights = TRUE, ev = NA, cng_dir, th_const = 1, noise_iter_max = 2, nob_init_min = 5, rmse = TRUE, use_last = TRUE, as_list = FALSE) {

  dim(x) <- c(nr, ny)
  dim(nob) <- c(nc, ny)
  
  v <- hilandyn_int_cpp(x, nob, nb, nc, ny, nr, yrs, foc_ind, cell_weights, n_times_eBias_of_mad, ev, cng_dir, th_const, noise_iter_max, nob_init_min, rmse, use_last)
  
  if (as_list) {
    valn <- c("D_MAX_YR", "D_MAX_MG", "D_FRS_YR", "D_FRS_MG", "D_LST_YR", "D_LST_MG", "D_NUM",
              "G_MAX_YR", "G_MAX_MG", "G_FRS_YR", "G_FRS_MG", "G_LST_YR", "G_LST_MG", "G_NUM", "N_GAP", "N_NOISE")
    vec1n <- c("RMSE")                                                                # vector (bands)
    vec2n <- c("LEN", "CPT_ID", "CPT_FOC", "NOISE", "MED_REL_MAG")                    # vector (years)
    matn <- c("EST", "SLO", "MAG", "REL_MAG")
    
    out_nm <- c(valn,
                paste0(rep(vec1n, each = nb), "_", seq(1, nb)),
                paste0(rep(vec2n, each = ny), "_", yrs),
                paste0(rep(matn, each = nb*ny), "_", seq(1, nb), "_", rep(yrs, each = nb)))
    
    if (length(v) == 1) v <- rep(NA, length(out_nm))

    names(v) <- out_nm
    
    return(list(D_MAX_YR = v[1], D_MAX_MG = v[2], D_FRS_YR = v[3], D_FRS_MG = v[4], D_LST_YR = v[5], D_LST_MG = v[6], D_NUM = v[7],
                G_MAX_YR = v[8], G_MAX_MG = v[9], G_FRS_YR = v[10], G_FRS_MG = v[11], G_LST_YR = v[12], G_LST_MG = v[13], G_NUM = v[14], N_GAP = v[15], N_NOISE = v[16], 
                RMSE = v[grep("RMSE", names(v))], 
                LEN = v[grep("LEN", names(v))], CPT_ID = v[grep("CPT_ID", names(v))], CPT_FOC = v[grep("CPT_FOC", names(v))], NOISE = v[grep("^NOISE", names(v))],
                MED_REL_MAG = v[grep("MED_REL_MAG", names(v))], 
                EST = matrix(v[grep("EST", names(v))], nb, ny), 
                SLO = matrix(v[grep("SLO", names(v))], nb, ny), 
                MAG = matrix(v[grep("^MAG", names(v))], nb, ny),
                REL_MAG = matrix(v[grep("^REL_MAG", names(v))], nb, ny)))
  }
  else {
    return(v)
  }
}