#' High-dimensional detection of Landscape Dynamics engine
#'
#' This is the internal function called by \code{\link{hilandyn_map}}.
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
#'   \item{D_MAX_MD}{Median among bands using values of \code{D_MAX} (single value).}
#'   \item{D_FST_MD}{Median among bands using values of \code{D_FST} (single value).}
#'   \item{G_MAX_MD}{Median among bands using values of \code{G_MAX} (single value).}
#'   \item{D_MAX_YR}{Year corresponding to \code{D_MAX_MD} (single value).}
#'   \item{D_FST_YR}{Year corresponding to \code{D_FST_MD} (single value).}
#'   \item{G_MAX_YR}{Year corresponding to \code{G_MAX_MD} (single value).}
#'   \item{D_MAX_DR}{Duration of the disturbance corresponding to \code{D_MAX_MD} (single value).}
#'   \item{D_FST_DR}{Duration of the disturbance corresponding to \code{D_FST_MD} (single value).}
#'   \item{G_MAX_DR}{Duration of the greening corresponding to \code{G_MAX_MD} (single value).}
#'   \item{D_MAX_ID}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_MAX_MD} (single value).}
#'   \item{D_FST_ID}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_FST_MD} (single value).}
#'   \item{G_MAX_ID}{Type of change (\code{CPT_ID}) of the greening corresponding to \code{G_MAX_MD} (single value).}
#'   \item{N_GAP}{Number of gaps in the time series, if any (single value).}
#'   \item{N_NOISE}{Number of years containing impulsive noise, if any (single value).}
#'   \item{RMSE}{Root mean square error of each band in the focal cell (one value per band).}
#'   \item{LEN}{Length of segments (one value per year).}
#'   \item{CPT_ID}{Type of change (one value per year). One of the following values: 101 (disturbance); 102 (greening); 9 (other change). Otherwise \code{NA} for no change.}
#'   \item{CPT_FOC}{Number of changepoints detected in the focal cell. The minimum value is zero and the maximum corresponds to the number of bands.}
#'   \item{NOISE}{Impulsive noise (one value per year).}
#'   \item{D_DUR}{Duration of disturbance (one value per year).}
#'   \item{G_DUR}{Duration of greening (one value per year).}
#'   \item{EST}{Estimated values of the bands in the focal cell (one value per year and band).}
#'   \item{SLO}{Slope of the linear segments of the bands in the focal cell (one value per year and band).}
#'   \item{MAG}{Magnitude in absolute terms of the bands in the focal cell (one value per year and band).}
#'   \item{MAG_REL}{Magnitude in relative terms of the bands in the focal cell (one value per year and band).}
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
#'                     cng_dir = cng_dir,
#'                     as_list = TRUE)
#'                 
#' # Plot results
#' dim(v) <- c(27, length(years))
#' par(mfrow=c(length(bands), 1))
#' 
#' for(i in seq_along(bands)) {
#'   plot(v[foc_ind[i],], type = "l", col = 2, lwd = 2, ylab = "", xlab = "", xaxt = "n")
#'   lines(out$EST[i,], type = "l", col = 4, lwd = 2, lty = 1)
#'   abline(v = which(out$CPT_ID %in% c(101, 201)), lty = 2, col = 2, lwd = 2)
#'   axis(1, at = seq_along(years), labels = years)
#'   title(main = paste(bands[i], "(RMSE =", out$RMSE[i,], ")"), adj = 0)
#' }
#'
#' @export

hilandyn_int <- function(x, nob, nb, nc, ny, nr, yrs, foc_ind, cell_weights = TRUE, ev = NA, cng_dir, th_const = 1, noise_iter_max = 2, nob_init_min = 5, rmse = TRUE, use_last = TRUE, as_list = FALSE) {
  
  dim(x) <- c(nr, ny)

  if (any(rowAlls(x, value = NA))) {
    return(ev)
  }
  
  dim(nob) <- c(nc, ny)

  cc <- !colAnyNAs(x)
  ncc <- sum(cc)
  
  if (ncc < ny) {
    
    gl <- rle_cpp(cc)
    gl <- gl[gl[, 1] == 0L, 2]
    
    if (any(gl > 1L)) {
      return(ev)
    }
    else {
      x <- x[, cc, drop = FALSE]
      nob <- nob[, cc, drop = FALSE]
    }
  }
  
  sd <- sd_est_cpp(x, n.times.eBias.of.mad)
  
  if (any(sd == 0L)) {
    return(ev)
  }
  
  if (cell_weights) {
    wgts <- wgt_sam_cpp(x, nb, nc)
  }
  else {
    wgts <- rep(1L, nr)
  }

  ntpa <- iter <- 0L
  cpt_esc <- tpav <- integer(dim(x)[2])
  
  if (dim(x)[2] < 5L) {
    return(ev)
  }
  
  y <- hd_bts_cpt_cpp(x, sd, nb, th_const, wgts, foc_ind - 1)
  
  if (rmse) {
    rmse <- rmse_cpp(x, y$est, foc_ind)
  }
  else {
    rmse <- rep(NA, nb)
  }

  if (noise_iter_max > 0L) {
    
    while (iter < noise_iter_max) {
      
      iter <- iter + 1L
      tpa <- 0L

      # search for the presence of noise
      if (length(y$cpt) > 0L) {
        
        cpt_in_ind <- !y$cpt %in% which(cpt_esc == 1)
        cpt_in <- y$cpt[cpt_in_ind]
        cpt_esc[cpt_in] <- 1L
        
        chk_mat <- cpt_cnd_cpp(x, y$est, nob, cpt_in, nb, nc, nob_init_min)

        # search for impulsive noise
        if (nrow(chk_mat) > 0L) {
          
          proc_cpt_res <- proc_cpt_cpp(x, y$cptind[, cpt_in_ind, drop = FALSE], chk_mat, cpt_in, sd, n.times.eBias.of.mad, wgts, nc, nb, th_const)  
          tpa <- proc_cpt_res$tpa
          
          if (any(tpa > 0L)) {
            ntpa <- sum(ntpa, length(tpa))
            tpav[tpa] <- 1L
            x <- proc_cpt_res$x
            y <- proc_cpt_res$y
          }
          else {
            break
          }
        }
        else {
          break
        }
      }
      else {
        break
      }
    }
  }
  
  if (ncc < ny) {
    ec <- which(cc == FALSE)
    
    # update cpts if there were data gaps
    if (length(y$cpt) > 0L) {
      fc <- seq_len(ny)[cc]
      
      for (i in seq_along(y$cpt)) {
        y$cpt[i] <- fc[y$cpt[i]]
      }
    }
    
    # fill missing columns
    for (i in seq_along(ec)) {
      eci <- ec[i]
      
      # process gaps occurring at the beginning
      if (eci == 1L) {
        
        x <- cbind(rep(NA, nr), x)
        y$est <- cbind(rep(NA, nr), y$est)
        tpav <- c(0L, tpav)
        
        if (any(y$cpt == eci + 1L | y$cpt == eci + 2L)) {
          # use the following value
          x[, eci] <- y$est[, eci] <- y$est[, eci + 1L]
        }
        else {
          # use extrapolated values
          x[, eci] <- y$est[, eci] <- y$est[, eci + 2L] + 2 * (y$est[, eci + 1L] - y$est[, eci + 2L])
        }
      }
      # or at the end
      else if (eci == ny) {
        x <- cbind(x, rep(NA, nr))
        y$est <- cbind(y$est, rep(NA, nr))
        tpav <- c(tpav, 0L)
        
        if (any(y$cpt == eci - 1L | y$cpt == eci - 2L)) {
          # use the preceding value
          x[, eci] <- y$est[, eci] <- y$est[, eci - 1L]
        }
        else {
          # use extrapolated values
          x[, eci] <- y$est[, eci] <- y$est[, eci - 2L] + 2 * (y$est[, eci - 1L] - y$est[, eci - 2L])
        }
      }
      else {
        seq1 <- seq_len(eci - 1L)
        seq2 <- eci:dim(x)[2]
        x <- cbind(x[, seq1, drop = FALSE], rep(NA, nr), x[, seq2, drop = FALSE])
        y$est <- cbind(y$est[, seq1, drop = FALSE], rep(NA, nr), y$est[, seq2,drop = FALSE])
        tpav <- c(tpav[seq1], 0L, tpav[seq2])
        
        if (any(y$cpt == eci + 1L)) {
          
          if (eci > 2L) {
            # fill gap using values extrapolated from the preceding segment
            x[, eci] <- y$est[, eci] <- y$est[, eci - 2L] + 2 * (y$est[, eci - 1L] - y$est[, eci - 2L])
          }
          else {
            # use the preceding value
            x[, eci] <- y$est[, eci] <- y$est[, eci - 1L]
          }
        }
        else {
          x[, eci] <- y$est[, eci] <- colMeans2(rbind(y$est[, eci - 1L], y$est[, eci + 1L]))
        }
      }
    }
  }
  
  # compute the number of segments
  nseg <- length(y$cpt) + 1L
  
  # compute the number of gaps
  ngap <- ny - ncc
  
  # assign the year to the gap occurrence
  tpav <- tpav * yrs
  
  # compute the indices at which every segment starts
  seg_beg <- c(1, y$cpt)
  
  # compute length of segments
  len_seg <- rep(NA, ny)
  
  for (i in seq_along(seg_beg)) {
    
    if (length(seg_beg) == 1L) {
      len <- ny
    }
    else if (i == length(seg_beg)) {
      len <- ny - seg_beg[i] + 1L
    }
    else {
      len <- seg_beg[i + 1] - seg_beg[i]
    }
    
    len_seg[seg_beg[i]] <- len
  }
  
  # compute slope values and place them at the beginning of each segment
  slo_seg <- matrix(NA, nr, ny)
  
  for (i in seq_along(seg_beg)) {
    
    if (len_seg[seg_beg[i]] == 1L) {
      slo <- rep(0L, nr)
    }
    else {
      slo <- y$est[, seg_beg[i] + 1L] - y$est[, seg_beg[i]]
    }
    
    slo_seg[, seg_beg[i]] <- slo
  }
  
  # compute the magnitude of change and identify the type of change
  mag <- mag_rel <- matrix(NA, nr, ny)
  cpt_id <- cpt_foc <- rep(NA, ny)
  d_mag_yr <- d_fst_yr <- g_mag_yr <- NA
  d_mag <- d_fst_md <- g_mag <- NA
  
  d_dur <- g_dur <- rep(NA, ny)
  d_mag_dr <- g_mag_dr <- d_fst_dr <- d_mag_id <- d_fst_id <- g_mag_id <- NA
  
  if (nseg > 1L) {
    
    hnb <- floor(nb/2)
    hnc <- foc_ind[1]

    for (i in seq_along(y$cpt)) {
      
      cpti <- y$cpt[i]
      cptind_cell <- matrix(y$cptind[, i], nc, nb)
      cptind_sum <- rowSums2(cptind_cell)
      cpt_foc[cpti] <- cptind_sum[hnc]

      mag1_act <- matrix(x[, cpti - 1L] - x[, cpti], nc, nb)
      mag1_est <- matrix(y$est[, cpti - 1L] - y$est[, cpti], nc, nb)
      
      cng_dir_act <- sign(mag1_act)
      cng_dir_est <- sign(mag1_est)
      
      sd_mat <- matrix(sd, nc, nb)
      cng_cnd <- sum(abs(mag1_est[hnc,]) < sd_mat[hnc,])

      cnd_d <- rowSums2(cng_dir_act == cng_dir & cng_dir_est == cng_dir)    
      cnd_g <- rowSums2(cng_dir_act == -cng_dir & cng_dir_est == -cng_dir)

      if (cpti == ny && !use_last) {
        cpt_id[cpti] <- 9L
      }

      else {  
        if (cnd_d[hnc] >= hnb) {
          cpt_id[cpti] <- 101L    # disturbance
          d_dur[cpti] <- len_seg[cpti]
          mag[, cpti] <- mag1_est
          mag_rel[, cpti] <- mag1_est / y$est[, cpti - 1L] * 100
        }
        else if (cnd_g[hnc] >= hnb) {    
          cpt_id[cpti] <- 102L    # greening
          g_dur[cpti] <- len_seg[cpti]
          mag[, cpti] <- mag1_est
          mag_rel[, cpti] <- mag1_est / y$est[, cpti - 1L] * 100
        }
        else {
          cpt_id[cpti] <- 9L
        }
      }
    }
    
    # find the highest-magnitude disturbance and its duration
    mag_med <- colMedians(abs(mag_rel))
    
    if (any(cpt_id %in% c(101))) {
      d_mag <- max(mag_med[cpt_id %in% c(101)])
      d_mag_ci <- which(mag_med %in% d_mag)
      d_mag_yr <- yrs[d_mag_ci]
      d_mag_dr <- d_dur[d_mag_ci]
      d_mag_id <- cpt_id[d_mag_ci]
      
      # find the first disturbance occurred
      ind <- which.max(cpt_id %in% c(101))
      d_fst_yr <- yrs[ind]
      d_fst_md <- mag_med[ind]
      d_fst_dr <- d_dur[ind]
      d_fst_id <- cpt_id[ind]
    }
    
    # find the highest-magnitude greening and its duration
    if (any(cpt_id %in% c(102))) {
      g_mag <- max(mag_med[cpt_id %in% c(102)])
      g_mag_ci <- which(mag_med %in% g_mag)
      g_mag_yr <- yrs[g_mag_ci]
      g_mag_dr <- g_dur[g_mag_ci]
      g_mag_id <- cpt_id[g_mag_ci]
    }
  }

  if (as_list) {
    v <- list(D_MAX_MD = d_mag, D_FST_MD = d_fst_md, G_MAG_MD = g_mag, D_MAX_YR = d_mag_yr,                                      # single values
              D_FST_YR = d_fst_yr, G_MAX_YR = g_mag_yr, D_MAX_DR = d_mag_dr, D_FST_DR = d_fst_dr,                                # single values
              G_MAX_DR = g_mag_dr, D_MAX_ID = d_mag_id, D_FST_ID = d_fst_id, G_MAX_ID = g_mag_id, N_GAP = ngap, N_NOISE = ntpa,  # single values
              RMSE = rmse,                                                                                                       # vectors (bands)
              LEN = len_seg, CPT_ID = cpt_id, CPT_FOC = cpt_foc, NOISE = tpav, D_DUR = d_dur, G_DUR = g_dur,                     # vectors (years)
              EST = y$est[foc_ind, , drop = FALSE], SLO = slo_seg[foc_ind, , drop = FALSE],                                      # matrices
              MAG = mag[foc_ind, , drop = FALSE], MAG_REL = mag_rel[foc_ind, , drop = FALSE])                                    # matrices                                                # matrices
  }
  else {
    v <- c(d_mag, d_fst_md, g_mag, d_mag_yr, d_fst_yr,      # single values
           g_mag_yr, d_mag_dr, d_fst_dr, g_mag_dr,          # single values
           d_mag_id, d_fst_id, g_mag_id, ngap, ntpa,        # single values
           rmse,                                            # vectors (bands)
           len_seg, cpt_id, cpt_foc, tpav, d_dur, g_dur,    # vectors (years)
           y$est[foc_ind, ], slo_seg[foc_ind, ],            # matrices
           mag[foc_ind, ], mag_rel[foc_ind, ])              # matrices
  }
  
  return(v)
}