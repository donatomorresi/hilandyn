#' High-dimensional detection of Landscape Dynamics
#'
#' This is the main function of the package.
#'
#' \code{hilandyn_map()} produces maps relative to landscape dynamics through the segmentation of high-dimensional Landsat time series. 
#' These latter include information from the spatial and spectral domains and are analysed using the High-dimensional Trend Segmentation (HiTS) procedure proposed by \insertCite{maeng2019adaptive;textual}{hilandyn}.
#' The HiTS procedure aims to detect changepoints in a piecewise linear signal where their number and location are unknown. Changes can occur in the intercept, slope or both of linear trends.
#' High-dimensional time series can include single or multiple spectral bands/indices, hereafter referred to as bands.
#' Impulsive noise, i.e. outliers in the time series, are removed through an iterative procedure.
#' One-year gaps in the time series are filled using either linear interpolation or extrapolation.
#' If input data, i.e. \code{sr_data}, \code{si_data}, and \code{nob_data}, are provided as paths to folders, one multiband raster in \emph{.tif} format per year is required. Raster file names should include the year such that they can be ordered by time.
#' Input rasters are read as \code{SpatRaster} objects created by the \pkg{terra} package.
#' 
#' @param sr_data character string. Surface reflectance data. Either the name of a \code{SpatRaster} object or the folder where rasters are stored. If \code{NULL} only spectral indices are used. See details for more information.
#' @param si_data character string. Spectral indices data. Either the name of a \code{SpatRaster} object or the folder where rasters are stored. If \code{NULL} only reflectance bands are used. See details for more information.
#' @param nob_data character string. Number of clear observations available per-pixel for producing reflectance composites. Either the name of a \code{SpatRaster} object or the folder where rasters are stored. If \code{NULL} a dummy value of 999 is used. See details for more information.
#' @param out_path character. The path where output rasters (in \emph{.tif} format) will be saved. Folders are created recursively if they do not exist. If \code{NULL} the output is a \code{SpatRaster} object created by the \pkg{terra} package.
#' @param sr_ind numeric vector. Indices of the layers containing reflectance data in single-date input rasters. If \code{NULL} all layers are used. Ignored if \code{sr_data} is a \code{SpatRaster}.
#' @param si_ind numeric vector. Indices of the layers containing spectral indices data in single-date input rasters. If \code{NULL} all layers are used. Ignored if \code{si_data} is a \code{SpatRaster}.
#' @param nob_ind numeric vector. Index of the layer containing data relative to the number of observations in single-date input rasters. If \code{NULL} the first layer of the input raster is used. Ignored if \code{nob_data} is a \code{SpatRaster}.
#' @param years numeric vector. Time interval (years) covered by the input rasters.
#' @param win_side integer. Width (in cells) of the spatial kernel used to extract input data from rasters. Must be an odd number.
#' @param cell_weights logical. Enable or disable the use of cell-based weights within the spatial kernel.
#' @param cng_dir numeric vector. Direction (either 1 or -1) of change in each spectral variable associated to the occurrence of a forest disturbance. It is computed as the difference between pre-disturbance and post-disturbance values.
#' @param th_const numeric. Constant value controlling the sensitivity to deviations from linearity in the HiTS procedure. Typical values are comprised in the interval \eqn{[0.7, 1.3]}.
#' @param nob_init_min integer. The minimum number of clear observations available per time step in the first two time steps of the time series.
#' @param noise_iter_max integer. The maximum number of iterations allowed for removing impulsive noise with the noise filter. The noise filter is disabled when the value is 0.
#' @param rmse logical. Determines whether to compute the root mean square error for each band in the focal cell. Original and estimated values are processed before noise filtering and gap filling.
#' @param use_last logical. Determines whether changepoints detected at the last time point are ignored or not.
#' @param expand logical. Whether to expand the raster by adding virtual rows and columns outside of it. New cells are filled with \code{NA}.
#' @param roi_vec character. Path to a spatial vector file containing polygons, e.g. a shapefile, of the region of interest. It is read as a \code{SpatVector} object using the \pkg{terra} package, and is used for masking input raster data.
#' @param clip_input logical. Whether to clip the input rasters using the extent of the region of interest.
#' @param n_copy integer. Number of copies of the input raster required for generating the output. Increasing this value will cause a reduction in the size of raster chucks that used during computation. Recommended values ranges between 1 and 4. 
#' @param cores integer. Number of CPU cores employed for parallelising the analysis.
#'
#' @return If \code{out_path} is \code{NULL}, a \code{SpatRaster} containing the following layers.
#'   \item{D_MAX_MD}{Median among bands using values of \code{D_MAX} (single layer).}
#'   \item{D_FST_MD}{Median among bands using values of \code{D_FST} (single layer).}
#'   \item{G_MAX_MD}{Median among bands using values of \code{G_MAX} (single layer).}
#'   \item{D_MAX_YR}{Year corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{D_FST_YR}{Year corresponding to \code{D_FST_MD} (single layer).}
#'   \item{G_MAX_YR}{Year corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{D_MAX_DR}{Duration of the disturbance corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{D_FST_DR}{Duration of the disturbance corresponding to \code{D_FST_MD} (single layer).}
#'   \item{G_MAX_DR}{Duration of the greening corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{D_MAX_ID}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{D_FST_ID}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_FST_MD} (single layer).}
#'   \item{G_MAX_ID}{Type of change (\code{CPT_ID}) of the greening corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{N_GAP}{Number of gaps in the time series, if any (single layer).}
#'   \item{N_NOISE}{Number of years containing impulsive noise, if any (single layer).}
#'   \item{RMSE}{Root mean square error of each band in the focal cell (one layer per band).}
#'   \item{LEN}{Length of segments (one layer per year).}
#'   \item{CPT_ID}{Type of change (one layer per year). One of the following values: 101 (disturbance); 102 (greening); 9 (other change). Otherwise \code{NA} for no change.}
#'   \item{CPT_FOC}{Number of changepoints detected in the focal cell. The minimum value is zero and the maximum corresponds to the number of bands.}
#'   \item{NOISE}{Impulsive noise (one layer per year).}
#'   \item{D_DUR}{Duration of disturbance (one layer per year).}
#'   \item{G_DUR}{Duration of greening (one layer per year).}
#'   \item{EST}{Estimated values of the bands in the focal cell (one layer per year and band).}
#'   \item{SLO}{Slope of the linear segments of the bands in the focal cell (one layer per year and band).}
#'   \item{MAG}{Magnitude in absolute terms of the bands in the focal cell (one layer per year and band).}
#'   \item{MAG_REL}{Magnitude in relative terms of the bands in the focal cell (one layer per year and band).}
#'
#' @author Donato Morresi, \email{donato.morresi@@unito.it}
#'
#' @seealso \code{\link{hilandyn_int}}
#'
#' @references
#' \insertRef{maeng2019adaptive}{hilandyn}
#'
#' @examples
#' library(hilandyn)
#'
#' # Load raster data
#' data(lnd_sr)
#' data(lnd_si)
#' lnd_sr <- terra::rast(lnd_sr)
#' lnd_si <- terra::rast(lnd_si)
#'
#' # Process data
#' rsout <- hilandyn_map(sr_data = "lnd_sr",
#'                       si_data = "lnd_si",
#'                       nob_data = NULL,
#'                       years = 1985:2020,
#'                       cng_dir = c(-1, -1, 1, 1),
#'                       cores = 1)
#'
#' # Plot the maximum disturbance magnitude and the corresponding year
#' terra::plot(rsout[["D_MAX_MD"]])
#' terra::plot(rsout[["D_MAX_YR"]])
#' 
#'
#' @export
#' @import Rcpp
#' @import RcppArmadillo
#' @import matrixStats
#' @import terra
#' @import foreach
#' @import future
#' @import doFuture
#' @import progressr
#' @importFrom lubridate seconds_to_period
#' @importFrom Rdpack reprompt


hilandyn_map <- function(sr_data, si_data, nob_data, out_path = NULL, sr_ind = NULL, si_ind = NULL, nob_ind = NULL, years, win_side = 3, cell_weights = TRUE, cng_dir, th_const = 1, noise_iter_max = 2, nob_init_min = 5, rmse = TRUE, use_last = TRUE, expand = TRUE, roi_vec = NULL, clip_input = FALSE, n_copy = 4, cores = 1) { 

  if (is.null(sr_data) && is.null(si_data)) {
    stop("missing input data")
  }
  
  if (win_side %% 2 == 0) {
    stop("window sides must be odd")
  }
  
  message("\nLoading raster data ...")

  if (!is.null(sr_data)) {
    
    if (!is.character(sr_data)) {
      stop("sr_data must be a character variable")
    }
      
    if (dir.exists(sr_data)) {
      f <- list.files(sr_data, ".tif$", full.names = TRUE)
      sr_data <- lapply(f, function(x) rast(x, lyrs = sr_ind))  # add the index for the number of observations
      sr_data <- do.call(c, sr_data)
    }
    else {
      sr_data <- get(sr_data)
      
      if (!inherits(sr_data, "SpatRaster")) {
        stop("can't load reflectance data")
      }
    }
  }
  
  if (!is.null(si_data)) {
    
    if (!is.character(si_data)) {
      stop("si_data must be a character variable")
    }
    
    if (dir.exists(si_data)) {
      f <- list.files(si_data, ".tif$", full.names = TRUE)
      si_data <- lapply(f, function(x) rast(x, lyrs = si_ind))
      si_data <- do.call(c, si_data)
    }
    else {
      si_data <- get(si_data)
      
      if (!inherits(si_data, "SpatRaster")) {
        stop("can't load spectral indices data")
      }
    }
  }
  
  if (!is.null(nob_data)) {
    
    if (!is.character(nob_data)) {
      stop("nob_data must be a character variable")
    }
    
    if (dir.exists(nob_data)) {
      f <- list.files(nob_data, ".tif$", full.names = TRUE)
      nob_data <- lapply(f, function(x) rast(x, lyrs = nob_ind))
      nob_data <- do.call(c, nob_data)
    }
    else {
      nob_data <- get(nob_data)
      
      if (!inherits(nob_data, "SpatRaster")) {
        stop("can't load data relative to the number of observations")
      }
    }
  }
      
  if (!is.null(sr_data) & !is.null(si_data)) {
    sr_mi <- matrix(seq_len(nlyr(sr_data)), ncol = nlyr(sr_data) / length(years), byrow = TRUE)
    si_mi <- matrix(seq_len(nlyr(si_data)), ncol = nlyr(si_data) / length(years), byrow = TRUE)
    rs <- lapply(seq_along(years), function(i) c(sr_data[[sr_mi[i,]]], si_data[[si_mi[i,]]]))
    rs <- do.call(c, rs)
  }
  else if (!is.null(sr_data)) {
    rs <- sr_data
  }
  else {
    rs <- si_data
  }
  
  if (!is.null(roi_vec)) {
    roi <- vect(roi_vec)
    
    if (crs(roi) != crs(rs)) {
      roi <- project(roi, crs(rs))
    }
    
    if (clip_input) {
      rs <- mask(crop(rs, roi), roi)
      nob_data <- mask(crop(nob_data, roi), roi)
    }
    else {
      rs <- mask(rs, roi)
      nob_data <- mask(nob_data, roi)
    }
  }

  ny <- length(years)
  nb <- nlyr(rs) / ny
  
  if (nb %% 1 != 0L) {
    stop("missing raster data for certain years")
  }
  
  if (length(cng_dir) != nb) {
    stop("length of cng_dir does not match the number of bands")
  }
  
  # Parameters for the focal computation
  nc <- win_side * win_side
  hw <- (win_side - 1) / 2
  nr <- nb * nc
  cell_ind <- matrix(seq_len(nr), nrow = nc)
  foc_cell <- win_side + floor(win_side / 2) + 1L
  foc_ind <- cell_ind[foc_cell, ]
  cng_dir <- rep(cng_dir, each = nc)

  # Output layer names
  valn <- c("D_MAX_MD", "D_FST_MD", "G_MAX_MD", "D_MAX_YR",          # single value
            "D_FST_YR", "G_MAX_YR", "D_MAX_DR", "D_FST_DR",
            "G_MAX_DR", "D_MAX_ID", "D_FST_ID", "G_MAX_ID",
            "N_GAP", "N_NOISE")
  vec1n <- c("RMSE")                                                 # vector (bands)
  vec2n <- c("LEN", "CPT_ID", "CPT_FOC", "NOISE", "D_DUR", "G_DUR")  # vector (years)
  matn <- c("EST", "SLO", "MAG", "MAG_REL")                          # matrix
  
  # Compute number of output layers
  nout <- (length(valn) + length(vec1n) * nb + length(vec2n) * ny + length(matn) * nb * ny)
  
  # Create empty vector for invalid pixels
  ev <- rep(NA, nout)
  
  # Set output parameters
  if (is.null(out_path)) {
    filename = ""
    overwrite = TRUE
    wopt = list()
  }
  else {
    if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
    filename = file.path(out_path, paste0("HILANDYN_", nb, "_BND_", nc, "_CEL_", years[1], "_", years[ny], ".tif"))
    overwrite = TRUE
    wopt=list(filetype='GTiff', datatype='FLT4S', NAflag = -9999, gdal=c("COMPRESS=LZW"))
  }

  out_nm <- c(valn,
              paste0(rep(vec1n, each = nb), "_", seq(1, nb)),
              paste0(rep(vec2n, each = ny), "_", years),
              paste0(rep(matn, each = nb*ny), "_", seq(1, nb), "_", rep(years, each = nb)))
  
  if (expand) {
    out <- rast(rs, nlyrs = nout)
    set.names(out, out_nm)
    rs <- extend(rs, hw)
    
    if (!is.null(nob_data)) {
      nob_data <- extend(nob_data, hw)
    }
  }
  else {
    out <- rast(nrows = nrow(rs) - 2*hw, 
                ncols = ncol(rs) - 2*hw, 
                nlyrs = nout, 
                crs = crs(rs), 
                extent = ext(rs) - hw * res(rs)[1], 
                resolution = res(rs))
    set.names(out, out_nm)
  }
  
  n_cols <- ncol(rs)
  readStart(rs)
  on.exit(readStop(rs))
  
  if (!is.null(nob_data)) {
    readStart(nob_data)
    on.exit(readStop(nob_data))
  }
  
  b <- writeStart(out, filename, overwrite, n = n_copy, wopt = wopt)
  tic <- proc.time()
  
  for (i in seq_len(b$n)) {
    
    message("\nProcessing chunk #", i, " of ", b$n)
    
    str_row <- b$row[i]
    n_rows <- b$nrows[i] + 2*hw
    v <- focal_values_cpp(rs@cpp$readValues(str_row-1, n_rows, 0, n_cols), c(n_rows, n_cols, nlyr(rs)), win_side)
    
    if (!is.null(nob_data)) {
      nob <- focal_values_cpp(nob_data@cpp$readValues(str_row-1, n_rows, 0, n_cols), c(n_rows, n_cols, nlyr(nob_data)), win_side)
    }
    else {
      nob = rep(999L, nc*ny)
    }
    
    if (cores > 1) {
      plan(multisession, workers = cores)
      on.exit(plan(sequential))
    }
    else {
      plan(sequential)
    }
    
    # Setup progress bar
    handlers(global = TRUE)
    on.exit(handlers(global = FALSE))
    p <- progressor(along = seq_len(ncol(v)))
    
    m <- foreach(j = seq_len(ncol(v)), .combine = 'cbind') %dofuture% {
      
      if (!is.null(nob_data)) {
        p()
        hilandyn_int(v[,j], nob[,j], nb, nc, ny, nr, years, foc_ind, cell_weights, ev, cng_dir, th_const, noise_iter_max, nob_init_min, rmse, use_last)
      }
      else {
        p()
        hilandyn_int(v[,j], nob, nb, nc, ny, nr, years, foc_ind, cell_weights, ev, cng_dir, th_const, noise_iter_max, nob_init_min, rmse, use_last)
      }
    }
    writeValues(out, t(m), str_row, b$nrows[i])
  }
  
  toc <- proc.time()
  elapsed <- round(seconds_to_period((toc - tic)[3]))
  message("\nProcessing took ", elapsed)

  out <- writeStop(out)
  gc()
  tmpFiles(remove = TRUE)
  return(out)
}
