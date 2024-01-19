#include <RcppArmadillo.h>
#include "hits.h"
#include "hits_misc.h"

// [[Rcpp::export]]
arma::umat rle_cpp(arma::uvec x) {
  
  // Get length(x)
  arma::uword n = x.n_elem;
  
  // Loop through and record value, length for runs
  arma::uword rows = 1;
  
  for (arma::uword a = 1; a < n; a++) {
    if (x(a) != x(a - 1)) rows += 1;
  }
  
  arma::umat out(rows, 4);
  arma::uword presval = x(0);
  arma::uword prespos = 0;
  arma::uword index = -1;
  
  for (arma::uword a = 1; a < n; a++) {
    
    arma::uword x_a = x(a);
    
    if (x_a != presval) {
      index += 1;
      out(index, 0) = presval;
      out(index, 1) = a - prespos;
      out(index, 2) = prespos;
      out(index, 3) = a - 1;  
      presval = x_a;
      prespos = a;
    }
  }
  
  index += 1;
  out(index, 0) = presval;
  out(index, 1) = n - prespos; 
  out(index, 2) = prespos;
  out(index, 3) = n - 1;
  
  return out;
}


// [[Rcpp::export]]
double mad_cpp(arma::vec x) {
  
  const double constant = 1.4826;
  arma::vec med(x.n_elem, arma::fill::value(median(x)));
  return arma::median(arma::abs(x - med)) * constant;
}


// [[Rcpp::export]]
arma::vec sd_est_cpp(arma::mat x, arma::vec mad_bias) {

  arma::uword n = x.n_cols - 2;

  double bias = 1 + mad_bias[n - 1] / n;

  arma::mat sd = arma::diff(x, 2, 1);

  arma::vec mad(sd.n_rows);

  for (arma::uword i = 0; i < sd.n_rows; i++) {
    mad[i] = mad_cpp(sd.row(i).t()) / bias / std::sqrt(6);
  }
  return mad;
}


// [[Rcpp::export]]
arma::vec wgt_sam_cpp(arma::mat x, arma::uword nb, arma::uword nc) {

  arma::vec sam_vec(nc, arma::fill::zeros);
  arma::uword hnc = (nc - 1) / 2;

  arma::uvec ind_vec = arma::regspace<arma::uvec>(0, x.n_rows - 1);
  arma::umat ind_mat = arma::reshape(ind_vec, nc, nb);

  for (arma::uword j = 0; j < nc; j++) {

    if (j == hnc) {
      continue;
    }
    else {
      arma::vec ngb_sam_vec(nb);  // for storing SAM values of different bands

      for (arma::uword n = 0; n < nb; n++) {
        
        double s = arma::dot(x.row(ind_mat(hnc, n)), x.row(ind_mat(j, n))) / (arma::norm(x.row(ind_mat(hnc, n))) * arma::norm(x.row(ind_mat(j, n))));
        
        if (s < -1) {
          s = -1;
        }
        
        if (s > 1) {
          s = 1;
        }

        ngb_sam_vec[n] = std::acos(s);
      }

      sam_vec[j] = arma::sum(ngb_sam_vec);
    }
  }

  arma::vec cell_wgt = arma::ones<arma::vec>(nc) - (sam_vec / arma::sum(sam_vec));
  return arma::repmat(cell_wgt, nb, 1);
}


// [[Rcpp::export]]
arma::mat focal_values_cpp(arma::vec d, arma::uvec dim, arma::uword w) {
  
  const int nrow = dim(0);
  const int ncol = dim(1);
  const arma::uword ncel = nrow * ncol;
  const arma::uword nlyr = dim(2);
  
  const int hw = (w  - 1) / 2;
  const arma::uword wcel = w * w;
  const arma::uword out_nr = nrow - 2*hw;
  const arma::uword out_nc = ncol - 2*hw;  
  
  const arma::uword nrow_val = wcel * nlyr;
  const arma::uword ncol_val = out_nr * out_nc;
  
  arma::mat val(nrow_val, ncol_val);
  
  if ((w / 2 == 0)) {
    Rcpp::Rcerr << "weights matrix must have uneven sides";
    return(val);
  }
  
  // create vector with start indices for the input vector
  arma::uvec in_str = arma::regspace<arma::uvec>(0, ncel, d.n_elem-1);
  arma::uvec lyr_ind = arma::repelem(in_str, w * w, 1); 
  
  arma::uvec cel_log(ncol);
  cel_log.subvec(0, out_nc - 1) = arma::ones<arma::uvec>(out_nc);
  cel_log = arma::repmat(cel_log, out_nr, 1);
  arma::uvec cel_ind = arma::find(cel_log == 1);  
  cel_ind = arma::repelem(cel_ind, w * w * nlyr, 1);
  
  arma::uvec smp_log(ncol);
  smp_log.subvec(0, w - 1) = arma::ones<arma::uvec>(w);
  smp_log = arma::repmat(smp_log, w, 1);
  arma::uvec smp_ind = arma::find(smp_log == 1);
  
  smp_ind = arma::repmat(smp_ind, nlyr, 1);
  smp_ind += lyr_ind;
  smp_ind = arma::repmat(smp_ind, ncol_val, 1);
  smp_ind += cel_ind;
  
  arma::mat out_mat(nrow_val, ncol_val, arma::fill::none);
  arma::uvec mat_ind = arma::regspace<arma::uvec>(0, smp_ind.n_elem - 1);
  out_mat(mat_ind) = d(smp_ind);
  
  return out_mat;
}


// [[Rcpp::export]]
arma::vec rmse_cpp(arma::mat x, arma::mat y, arma::uvec foc_ind) {
  
  arma::vec rmse = arma::sqrt(arma::sum(arma::square(x.rows(foc_ind) - y.rows(foc_ind)), 1) / x.n_cols);
  return rmse;
}


// [[Rcpp::export]]
arma::umat cpt_cnd_cpp(arma::mat x, arma::mat y, arma::mat nob, arma::uvec cpt, arma::uword nb, arma::uword nc, arma::uword nob_init_min) {
  
  cpt -= 1;                          // for compatibility between C++ and R
  arma::uword ncol = 4 + x.n_cols;
  arma::umat chk_mat(0, ncol);       // information on time points (length, start, end, condition for verification, indices of the observations to be removed)
  
  if (cpt.n_elem > 0) {
    
    arma::uvec col_ind = arma::regspace<arma::uvec>(0, x.n_cols - 1);
    arma::uvec row_ind = arma::regspace<arma::uvec>(0, x.n_rows - 1);
    arma::umat ind_mat = arma::reshape(row_ind, nc, nb);
    
    arma::vec nob_vec = arma::median(nob, 0).t();   // median number of observations within the spatial kernel
    arma::uvec nob_init_cnd = nob_vec < nob_init_min;
    
    // partition into sequences of changepoints
    arma::uvec cpt_vec(x.n_cols);
    cpt_vec(cpt) = arma::ones<arma::uvec>(cpt.n_elem);
    arma::umat seq_mat = rle_cpp(cpt_vec);
    seq_mat.shed_rows(arma::find(seq_mat.col(0) == 0));
    
    // store information in a matrix
    chk_mat.insert_rows(0, seq_mat.n_rows);
    chk_mat.col(0) = seq_mat.col(1);
    chk_mat.col(1) = seq_mat.col(2);
    chk_mat.col(2) = seq_mat.col(3);
    
    // iterate over sequences
    for (arma::uword i = 0; i < chk_mat.n_rows; i++) {
      
      arma::uword ind_1 = chk_mat(i, 1);       // time index of the first cpt
      arma::uword ind_2 = chk_mat(i, 2);       // time index of the last cpt
      arma::uword ind_0 = ind_1 - 1;           // time index preceding the cpt
      arma::uvec col_seq_ind = col_ind.subvec(ind_0, ind_2);
      
      if (ind_1 == x.n_cols - 1) {
        continue;
      }
      
      // Test conditions based on the number of observations
      arma::uvec low_nob_init_seq_ind = arma::find(nob_init_cnd.subvec(ind_0, ind_2) == 1);

      if ((nob_init_cnd[0] == 1 || nob_init_cnd[1] == 1) && ind_1 == 1) {
        
        arma::uvec low_nob_init_ind = col_seq_ind(low_nob_init_seq_ind);
        arma::uvec row_i = {i};
        chk_mat(row_i, arma::regspace<arma::uvec>(0, low_nob_init_ind[low_nob_init_ind.n_elem - 1]) + 4).ones();
        chk_mat(i, 3) = 1;
        continue;
      }
      else {
        
        arma::mat dst_mat(nb, col_seq_ind.n_elem);

        for (arma::uword j = 0; j < col_seq_ind.n_elem; j++) {
          
          if (col_seq_ind[j] == 0 || col_seq_ind[j] == x.n_cols - 1) {
            continue;
          }
          
          // Euclidean distance from reference values
          if (j == 0) {     // use original values
            
            for (arma::uword k = 0; k < nb; k++) {
              
              arma::uvec mat_ind_1 = ind_mat.col(k) + (x.n_rows * col_seq_ind[j]);
              arma::uvec dst_ref_ind;
              
              if (ind_0 == 0) {
                dst_ref_ind = ind_mat.col(k) + (x.n_rows * (ind_2 + 1));
              }
              else {
                dst_ref_ind = ind_mat.col(k) + (x.n_rows * (ind_0 - 1));
              }
              
              dst_mat(k, j) = arma::norm(y(dst_ref_ind) - x(mat_ind_1), 2);
            }
          }
          else {
            for (arma::uword k = 0; k < nb; k++) {
              
              arma::uvec mat_ind_1 = ind_mat.col(k) + (x.n_rows * col_seq_ind[j]);
              arma::uvec dst_ref_ind;

              if (ind_0 == 0) {
                dst_ref_ind = ind_mat.col(k) + (x.n_rows * (ind_2 + 1));
              }
              else {
                dst_ref_ind = ind_mat.col(k) + (x.n_rows * (ind_0 - 1));
              }

              dst_mat(k, j) = arma::norm(y(dst_ref_ind) - y(mat_ind_1), 2);
            }
          }
        }
        
        if (all(dst_mat.col(1) > dst_mat.col(0))) {
          dst_mat.shed_col(0);
          col_seq_ind.shed_row(0);
        }
        
        arma::uvec dst_ind_max = arma::unique(arma::index_max(dst_mat, 1));
        arma::uvec pnt_dst_max = col_seq_ind(dst_ind_max);
        arma::uvec pnt_dst_max_cnd(pnt_dst_max.n_elem);
        
        for (arma::uword l = 0; l < pnt_dst_max.n_elem; l++) {
          
          arma::uword i0 = {pnt_dst_max[l] - 1};
          arma::uword i1 = {pnt_dst_max[l]};
          arma::uword i2 = {pnt_dst_max[l] + 1};

          arma::uvec sgn(nb);
          arma::vec ed1(nb);
          arma::vec ed2(nb);

          if (i1 == ind_0) {
            for (arma::uword k = 0; k < nb; k++) {
              
              arma::uvec mat_ind_0 = ind_mat.col(k) + x.n_rows * i0;
              arma::uvec mat_ind_1 = ind_mat.col(k) + x.n_rows * i1;
              arma::uvec mat_ind_2 = ind_mat.col(k) + x.n_rows * i2;
              
              sgn[k] = arma::sum(arma::sign(x(mat_ind_0) - x(mat_ind_1)) != arma::sign(x(mat_ind_1) - x(mat_ind_2)));
              ed1[k] = arma::norm(x(mat_ind_1) - x(mat_ind_2));
              ed2[k] = arma::norm(x(mat_ind_0) - x(mat_ind_2));
            }
          }
          else {
            for (arma::uword k = 0; k < nb; k++) {
              
              arma::uvec mat_ind_0 = ind_mat.col(k) + x.n_rows * i0;
              arma::uvec mat_ind_1 = ind_mat.col(k) + x.n_rows * i1;
              arma::uvec mat_ind_2 = ind_mat.col(k) + x.n_rows * i2;
              
              sgn[k] = arma::sum(arma::sign(y(mat_ind_0) - y(mat_ind_1)) != arma::sign(y(mat_ind_1) - y(mat_ind_2)));
              ed1[k] = arma::norm(y(mat_ind_1) - y(mat_ind_2));
              ed2[k] = arma::norm(y(mat_ind_0) - y(mat_ind_2));
            }
          }
          arma::uvec cnd = ed2 < ed1 && sgn > 0;
          pnt_dst_max_cnd[l] = arma::sum(cnd) == 0;
        }
        
        pnt_dst_max.shed_rows(arma::find(pnt_dst_max_cnd == 1));
        
        arma::uvec ind_i = { i };
        chk_mat(ind_i, pnt_dst_max + 4).ones();
      }
    }
    chk_mat.shed_rows(arma::find(arma::sum(chk_mat.cols(4, chk_mat.n_cols - 1), 1) == 0));
  }
  return chk_mat;
}


// [[Rcpp::export]]
Rcpp::List proc_cpt_cpp(arma::mat x, arma::umat cptind, arma::umat chk_mat, arma::uvec cpt, arma::vec sd, arma::vec mad_bias, arma::vec wgts, arma::uword nc, arma::uword nb, double th_const) {
  
  cpt -= 1;
  arma::uvec tpa(x.n_cols);
  
  arma::umat row_ind_mat(nc, nb);
  arma::uvec ind_all = arma::regspace<arma::uvec>(0, x.n_rows - 1);
  row_ind_mat(ind_all) = ind_all;
  
  arma::uword hnc = (nc - 1) / 2;
  arma::uvec foc_ind = arma::trans(row_ind_mat.row(hnc));
  
  arma::umat cptind_mat(x.n_rows, x.n_cols);
  cptind_mat.cols(cpt) = cptind;
  
  for (arma::uword i = 0; i < chk_mat.n_rows; i++) {
    
    if (chk_mat(i, 3) == 1) {   // skip sequences that cannot be verified (at the beginning)
      continue;
    }
    
    arma::urowvec cnd_vec = chk_mat(i, arma::span(4, chk_mat.n_cols - 1));
    arma::uvec col_ind = arma::find(cnd_vec == 0);     // time indices to be retained for testing
    arma::uvec rm_ind = arma::find(cnd_vec == 1);      // time indices to be removed
    
    arma::umat cptind_i(nc, nb);
    arma::uvec cptind_vec = cptind_mat.col(chk_mat(i, 1));
    
    if (all(cptind_vec == 0)) {
      cptind_i(arma::regspace<arma::uvec>(0, cptind_mat.n_rows - 1)) = arma::ones<arma::uvec>(cptind_mat.n_rows);
    }
    else {
      cptind_i(arma::regspace<arma::uvec>(0, cptind_mat.n_rows - 1)) = cptind_vec;
    }
    
    // band selection approach
    arma::vec cell_cnt = arma::conv_to<arma::vec>::from(arma::sum(cptind_i, 0));
    arma::uvec band_ind = arma::find(cell_cnt >= arma::as_scalar(arma::median(cell_cnt)));
    arma::uvec row_ind = arma::vectorise(row_ind_mat.cols(band_ind));

    arma::uvec foc_ind_new(band_ind.n_elem);
    arma::umat row_ind_mat_new(nc, band_ind.n_elem);
    arma::uvec ind_all_new = arma::regspace<arma::uvec>(0, row_ind_mat_new.n_elem - 1);
    row_ind_mat_new(ind_all_new) = ind_all_new;
    foc_ind_new = arma::vectorise(row_ind_mat_new.row(hnc));

    arma::vec sd_new = sd_est_cpp(x(row_ind, col_ind), mad_bias);

    if (any(sd_new) == 0) {
      sd_new = sd(row_ind);
    }

    Rcpp::List y = hd_bts_cpt_cpp(x(row_ind, col_ind), sd_new, band_ind.n_elem, th_const, wgts(row_ind), foc_ind_new);

    arma::uword n_pre;
    if (rm_ind[0] < chk_mat(i, 1)) {
      n_pre = chk_mat(i, 1) - rm_ind[0];
    }
    else {
      n_pre = 0;
    }
    
    arma::uvec cpt_chk = { chk_mat(i, 1) };
    cpt_chk -= n_pre;
    arma::uvec cpt_tst = Rcpp::as<arma::uvec>(y["cpt"]) - 1;           // for compatibility between R and C++

    arma::vec match_cpt_ind = match_cpp(cpt_tst, cpt_chk);
    match_cpt_ind = match_cpt_ind(arma::find_finite(match_cpt_ind));

    if (match_cpt_ind.n_elem == 0) {
      chk_mat(i, 3) = 1;
    }
  }
  
  chk_mat.shed_rows(arma::find(chk_mat.col(3) == 0));
  Rcpp::List y(3);
  
  if (chk_mat.n_rows > 0) {
    
    for (arma::uword j = 0; j < chk_mat.n_rows; j++) {
      
      arma::uvec rm_ind = arma::find(chk_mat(j, arma::span(4, chk_mat.n_cols - 1)) == 1);

      if (rm_ind[0] == 0) {
        x.cols(rm_ind) = arma::repmat(arma::mean(x.cols(rm_ind[rm_ind.n_elem - 1] + 1, rm_ind[rm_ind.n_elem - 1] + 2), 1), 1, rm_ind.n_elem);
      }
      else if (rm_ind[rm_ind.n_elem - 1] == x.n_cols - 1) {     // extrapolate at the end
        x.cols(rm_ind) = arma::repmat(arma::mean(x.cols(rm_ind[0] - 1, rm_ind[0] - 2), 1), 1, rm_ind.n_elem);
      }
      else {                    // interpolate
        x.cols(rm_ind) = arma::repmat(arma::mean(arma::join_horiz(x.col(rm_ind[0] - 1), x.col(rm_ind[rm_ind.n_elem - 1] + 1)), 1), 1, rm_ind.n_elem);
      }
      
      tpa(rm_ind).ones();
    }
    
    y = hd_bts_cpt_cpp(x, sd, nb, th_const, wgts, foc_ind);
  }
  
  tpa = arma::find(tpa == 1);
  tpa += 1;
  
  return Rcpp::List::create(Rcpp::Named("x") = x,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("tpa") = tpa);
}
