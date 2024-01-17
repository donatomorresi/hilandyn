#' High-dimensional detection of Landscape Dynamics
#'
#' This is the main function of the package.
#'
#' \code{hilandyn_map()} produces maps relative to landscape dynamics through the segmentation of high-dimensional Landsat time series. 
#' These latter include information in the spatial and spectral dimensions and are analysed using the High-dimensional Trend Segmentation (HiTS) procedure proposed by \insertCite{maeng2019adaptive;textual}{hilandyn}.
#' The HiTS procedure aims to detect changepoints in a piecewise linear signal where their number and location are unknown. Changes can occur in the intercept, slope or both of linear trends.
#' Time series can include single or multiple spectral bands/indices, hereafter referred to as bands.
#' Impulsive noise, i.e. outliers in the time series, are removed through an iterative procedure.
#' One-year gaps in the time series are filled using either linear interpolation or extrapolation.
#' If input data are provided as folders, one multiband raster in \emph{.tif} format per year should be present within each folder. Moreover, filenames should contain the year such that rasters can be ordered by year.
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
#' @param use_last logical. Determines whether changepoints detected at the last time point are ignored or not.
#' @param expand logical. Whether to expand the raster by adding virtual rows and columns outside of it. New cells are filled with \code{NA}.
#' @param roi_vec character. Path to a spatial vector file containing polygons, e.g. a shapefile, of the region of interest. It is read as a \code{SpatVector} object using the \pkg{terra} package, and is used for masking input raster data.
#' @param clip_input logical. Whether to clip the input rasters using the extent of the region of interest.
#' @param n_copy integer. Number of copies of the input raster required for generating the output. Increasing this value will cause a reduction in the size of raster chucks that used during computation. Recommended values ranges between 1 and 4. 
#' @param cores integer. Number of CPU cores employed for parallelising the analysis.
#'
#' @return If \code{out_path} is \code{NULL}, a \code{SpatRaster} containing the following layers.
#'   \item{EST}{Estimated values (one layer for each year and band).}
#'   \item{CPT}{Detected changepoints (one layer for each year and band).}
#'   \item{SLO}{Slope of the linear segments (one layer for each year and band).}
#'   \item{MAG}{Magnitude in absolute terms (one layer for each year and band).}
#'   \item{MAG_REL}{Magnitude in relative terms (one layer for each year and band).}
#'   \item{LEN}{Length of segments (one layer per year).}
#'   \item{CPT_ID}{Type of change (one layer per year). One of the following values: 101 (abrupt disturbance); 102 (abrupt greening); 201 (gradual disturbance); 202 (gradual greening); 9 (other change). Otherwise \code{NA} for no change.}
#'   \item{CPT_FOC}{Number of changepoints detected by the HiTS procedure in the focal cell. The minimum value is zero and the maximum corresponds to the number of bands.}
#'   \item{NOISE}{Impulsive noise (one layer per year).}
#'   \item{D_DUR}{Duration of disturbance (one layer per year).}
#'   \item{G_DUR}{Duration of greening (one layer per year).}
#'   \item{D_MAX}{Maximum disturbance change magnitude (in relative terms) throughout the time series (one layer per band).}
#'   \item{D_FST}{Change magnitude (in relative terms) associated with the first disturbance detected within the time series (one layer per band).}
#'   \item{G_MAX}{Maximum greening change magnitude (in relative tems) throughout the time series (one layer per band).}
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
#'
#' @author Donato Morresi, \email{donato.morresi@@gmail.com}
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
#' rsout <- hilandyn(sr_data = "lnd_sr",
#'                   si_data = "lnd_si",
#'                   years = 1985:2020,
#'                   cng_dir = c(-1, -1, 1, 1),
#'                   cores = parallel::detectCores() - 1)
#'
#' # Plot the maximum disturbance magnitude
#' terra::plot(rsout[["D_MAX_MD"]])
#'
#' @export
#' @import Rcpp
#' @import RcppArmadillo
#' @import matrixStats
#' @import terra
#' @import foreach
#' @import future
#' @import doFuture
#' @import cli
#' @import progressr
#' @importFrom Rdpack reprompt


hilandyn_map <- function(sr_data, si_data, nob_data, out_path = NULL, sr_ind = NULL, si_ind = NULL, nob_ind = NULL, years, win_side = 3, cell_weights = TRUE, cng_dir, th_const = 1, noise_iter_max = 2, nob_init_min = 5, use_last = TRUE, expand = TRUE, roi_vec = NULL, clip_input = FALSE, n_copy = 4, cores = 1) { 

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
  valn <- c("D_MAX_MD", "D_FST_MD", "G_MAX_MD", "D_MAX_YR",          # single values
            "D_FST_YR", "G_MAX_YR", "D_MAX_DR", "D_FST_DR",
            "G_MAX_DR", "D_MAX_ID", "D_FST_ID", "G_MAX_ID",
            "N_GAP", "N_NOISE")
  vecn <- c("LEN", "CPT_ID", "CPT_FOC", "NOISE", "D_DUR", "G_DUR")   # long vectors (years)
  matn <- c("EST", "SLO", "MAG", "MAG_REL")                          # matrices
  
  # Compute number of output layers
  nout <- (length(valn) + length(vecn) * ny + length(matn) * nb * ny)
  
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

  out_nm <- c(valn, paste0(rep(vecn, each = ny), "_", years), 
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
  
  # Set rasters properties
  n_cols <- ncol(rs)

  # Begin reads
  readStart(rs)
  on.exit(readStop(rs))
  
  if (!is.null(nob_data)) {
    readStart(nob_data)
    on.exit(readStop(nob_data))
  }
  
  b <- writeStart(out, filename, overwrite, n = n_copy, wopt = wopt)
  
  for (i in seq_len(b$n)) {
    
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
    
    p <- progressor(along = ncol(v))
     
    m <- foreach(j = seq_len(ncol(v)), .combine = 'cbind') %dofuture% {
    
      if (!is.null(nob_data)) {
        hilandyn_int(v[,j], nob[,j], nr, ny, nb, nc, years, foc_ind, cell_weights, ev, cng_dir, th_const, noise_iter_max, nob_init_min, use_last)
      }
      else {
        hilandyn_int(v[,j], nob, nr, ny, nb, nc, years, foc_ind, cell_weights, ev, cng_dir, th_const, noise_iter_max, nob_init_min, use_last)
      }
      p(sprintf("Chunk ", j))
    }
    writeValues(out, t(m), str_row, b$nrows[i])
  }

  out <- writeStop(out)
  gc()
  tmpFiles(remove = TRUE)
  return(out)
}
