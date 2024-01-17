#############################################################################################################
## Find potential point anomalies
cpt_cnd <- function(cpt, x, sd, dir, nb, nc) {
  
  cpt_chk <- numeric(dim(x)[2])
  cpt <- cpt[cpt != dim(x)[2]]
  
  if (length(cpt) > 0) {
    
    hrows <- round(dim(x)[1] / 2)
    
    # consecutive changepoints
    seq_cpt <- cpt[diff2(cpt) == 1]
    
    if (length(seq_cpt) > 0) {
      
      # cpt at the beginning of the sequence
      cnd <- sign(x[, seq_cpt[1] - 1L] - x[, seq_cpt[1]]) != dir
      cnd_sum<- rowSums(matrix(cnd, ncol = nb))
      
      if (any(cnd_sum > 0)) {
        cpt_chk[seq_cpt[1]] <- 1L
      }
    }
    
    # changepoints offsets
    cpt_off <- unique(c(cpt, cpt - 1L))
    cpt_off <- cpt_off[!(cpt_off %in% c(1, which(cpt_chk == 1)))]
    
    for (i in seq_along(cpt_off)) {
      
      cpti <- cpt_off[i]
      
      v0 <- x[, cpti - 1]
      v1 <- x[, cpti]
      v2 <- x[, cpti + 1]
      
      dst0 <- acos(sum(v0 * v2) / (sqrt(sum(v0 ^ 2)) * sqrt(sum(v2 ^ 2))))
      dst1 <- acos(sum(v1 * v2) / (sqrt(sum(v1 ^ 2)) * sqrt(sum(v2 ^ 2))))
      
      dst2 <- sqrt(sum(abs(v0 / sd - v2 / sd) ^ 2))
      dst3 <- sqrt(sum(abs(v1 / sd - v2 / sd) ^ 2))
      
      mag1 <- x[, cpti - 1] - x[, cpti]
      mag2 <- x[, cpti] - x[, cpti + 1]
      
      cnd_sum <- sum(sign(mag1) != sign(mag2))
      
      if (cnd_sum > hrows && (dst0 < dst1 || dst2 < dst3)) {
        cpt_chk[cpti] = 1;
      }
    }
  }
  return(which(cpt_chk == 1)) 
  
}


##############################################################################################################
## Process potential point anomalies
proc_cpt <- function(x, mad_bias, cpt, th_const, wgts, dir, nb) {
  
  tpa <- numeric(length(cpt))
  
  for (i in seq_along(cpt)) {
    
    sd <- sd_est_cpp(x[, -cpt[i], drop = FALSE], mad_bias)
    yp <- hd_bts_cpt_cpp(x[, -cpt[i], drop = FALSE], sd, nb, th_const, wgts)
    
    cpt_off <- c(cpt[i], cpt[i] - 1, cpt[i] + 1)
    cpt_tst <- yp$cpt
    
    if (any(cpt_tst %in% cpt_off)) {
      next
    }
    else {
      tpa[i] <- cpt[i]
    }
  }
    
  return(tpa[tpa != 0])
}












