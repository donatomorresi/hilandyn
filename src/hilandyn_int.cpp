#include <RcppArmadillo.h>
#include "hits.h"
#include "hits_misc.h"
#include "hilandyn_misc.h"

// [[Rcpp::export]]
arma::vec hilandyn_int_cpp(arma::mat x, arma::mat nob, arma::uword nb, arma::uword nc, arma::uword ny, arma::uword nr, arma::uvec yrs, arma::uvec foc_ind, arma::uword cell_weights, arma::vec n_times_eBias_of_mad, arma::vec ev, arma::ivec cng_dir, double th_const, arma::uword noise_iter_max, arma::uword nob_init_min, arma::uword rmse, arma::uword use_last) {
  
  x.reshape(nr, ny);
  nob.reshape(nc, ny);
  
  arma::umat nfnt_x(nr, ny, arma::fill::zeros);
  nfnt_x(arma::find_nonfinite(x)).ones();
  
  if (any(arma::sum(nfnt_x, 1) == ny)) {
    return ev;
  }

  arma::uvec empty_col_cnt = arma::sum(nfnt_x, 0).t();
  arma::uvec cc_vec(ny, arma::fill::ones);
  cc_vec(arma::find(empty_col_cnt > 0)).zeros();
  arma::uword ncc = arma::sum(cc_vec);
  
  if (ncc < ny) {
    arma::umat gl_mat = rle_cpp(cc_vec);
    gl_mat.shed_rows(arma::find(gl_mat.col(0) == 1));
    arma::uvec gl_vec = gl_mat.col(1);
    
    if (any(gl_vec > 1)) {
      return ev;
    }
    else {
      x.shed_cols(arma::find(cc_vec == 0));
      nob.shed_cols(arma::find(cc_vec == 0));
    }
  }
  
  if (x.n_cols < 6) {
    return ev;
  }
  
  arma::vec sd = sd_est_cpp(x, n_times_eBias_of_mad);
  
  if (any(sd == 0)) {
    return ev;
  }
  
  arma::vec wgts = arma::ones<arma::vec>(nr);
  
  if (cell_weights == 1) {
    wgts = wgt_sam_cpp(x, nb, nc);
  }
  
  arma::uword ntpa = 0, iter = 0;
  arma::uvec cpt_esc(x.n_cols);
  arma::vec tpav(x.n_cols, arma::fill::value(arma::datum::nan));
  
  Rcpp::List y = hd_bts_cpt_cpp(x, sd, nb, th_const, wgts, foc_ind);
  arma::uvec cpt = Rcpp::as<arma::uvec>(y["cpt"]);
  arma::umat cptind = Rcpp::as<arma::umat>(y["cptind"]);
  cpt -= 1;
  
  if (noise_iter_max > 0) {
    while (iter < noise_iter_max) {
      iter += 1;

      // search for the presence of noise
      if (cpt.n_elem > 0) {
        arma::vec cpt_in_ind0 = match_cpp(cpt, arma::find(cpt_esc == 1));
        arma::uvec cpt_in_ind = arma::find_nonfinite(cpt_in_ind0);
        arma::uvec cpt_in = cpt(cpt_in_ind);
        cpt_esc(cpt_in).ones();
        arma::umat chk_mat = cpt_cnd_cpp(x, Rcpp::as<arma::mat>(y["est"]), nob, cpt_in, nb, nc, nob_init_min);
        
        // search for impulsive noise
        if (chk_mat.n_rows > 0) {
          x = proc_cpt_cpp(x, cptind.cols(cpt_in_ind), chk_mat, cpt_in, sd, n_times_eBias_of_mad, wgts, nc, nb, th_const);
          arma::rowvec tpa = x.row(nr);
          arma::uvec tpa_ind = arma::find(tpa == 1);
          x.shed_row(nr);
          
          if (tpa_ind.n_elem > 0) {
            ntpa = ntpa + tpa_ind.n_elem;
            tpav(tpa_ind).ones();
            y = hd_bts_cpt_cpp(x, sd, nb, th_const, wgts, foc_ind);
            cpt = Rcpp::as<arma::uvec>(y["cpt"]);
            cptind = Rcpp::as<arma::umat>(y["cptind"]);
            cpt -= 1;
          }
          else {
            break;
          }
        }
        else {
          break;
        }
      }
      else {
        break;
      }
    }
  }

  arma::mat y_est = Rcpp::as<arma::mat>(y["est"]);

  if (ncc < ny) {
    arma::uvec ec_ind = arma::find(cc_vec == 0); 
    arma::uvec cc_ind = arma::find(cc_vec == 1);
    arma::uword ncp = cpt.n_elem;

    // update cpts if there were data gaps
    if (ncp > 0) {
      arma::uvec fc = arma::regspace<arma::uvec>(0, ny - 1);
      fc.shed_rows(ec_ind);
      
      for (arma::uword i = 0; i < ncp; i++) {
        cpt(i) = fc(cpt(i));
      }
    }
    
    // fill missing columns
    arma::vec two_vec(nr, arma::fill::value(2));
    
    for (arma::uword i = 0; i < ec_ind.n_elem; i++) {
      arma::uword eci = ec_ind[i];
      
      // process gaps occurring at the beginning
      if (eci == 0) {
        x = arma::join_horiz(arma::zeros<arma::vec>(nr), x);
        y_est = arma::join_horiz(arma::zeros<arma::vec>(nr), y_est);
        tpav = arma::join_vert(arma::vec(1, arma::fill::value(arma::datum::nan)), tpav);  
        
        if (any(cpt == eci + 1 || cpt == eci + 2)) {
          // use the following value
          x.col(eci) = y_est.col(eci + 1);
          y_est.col(eci) = y_est.col(eci + 1);
        }
        else {
          // use extrapolated values
          x.col(eci) = y_est.col(eci + 2) + two_vec % (y_est.col(eci + 1) - y_est.col(eci + 2));
          y_est.col(eci) = y_est.col(eci + 2) + two_vec % (y_est.col(eci + 1) - y_est.col(eci + 2));
        }
      }
      // or at the end
      else if (eci == ny - 1) {
        x = arma::join_horiz(x, arma::zeros<arma::vec>(nr));
        y_est = arma::join_horiz(y_est, arma::zeros<arma::vec>(nr));
        tpav = arma::join_vert(tpav, arma::vec(1, arma::fill::value(arma::datum::nan)));
        
        if (any(cpt == eci - 1 || cpt == eci - 2)) {
          // use the preceding value
          x.col(eci) = y_est.col(eci - 1);
          y_est.col(eci) = y_est.col(eci - 1);
        }
        else {
          // use extrapolated values
          x.col(eci) = y_est.col(eci - 2) + two_vec % (y_est.col(eci - 1) - y_est.col(eci - 2));
          y_est.col(eci) = y_est.col(eci - 2) + two_vec % (y_est.col(eci - 1) - y_est.col(eci - 2));
        }
      }
      else {
        arma::uvec seq1 = arma::regspace<arma::uvec>(0, eci - 1);
        arma::uvec seq2 = arma::regspace<arma::uvec>(eci, x.n_cols - 1);
        x = arma::join_horiz(x.cols(seq1), arma::zeros<arma::vec>(nr), x.cols(seq2));
        y_est = arma::join_horiz(y_est.cols(seq1), arma::zeros<arma::vec>(nr), y_est.cols(seq2));
        tpav = arma::join_vert(tpav(seq1), arma::vec(1, arma::fill::value(arma::datum::nan)), tpav(seq2));
        
        if (any(cpt == eci + 1)) {
          if (eci > 1) {
            // fill gap using values extrapolated from the preceding segment
            x.col(eci) = x.col(eci - 2) + two_vec % (x.col(eci - 1) - x.col(eci - 2));
            y_est.col(eci) = y_est.col(eci - 2) + two_vec % (y_est.col(eci - 1) - y_est.col(eci - 2));
          }
          else {
            // use the preceding value
            x.col(eci) = y_est.col(eci - 1);
            y_est.col(eci) = y_est.col(eci - 1);
          }
        }
        else {
          arma::mat ngb_cols = arma::join_horiz(y_est.col(eci - 1), y_est.col(eci + 1));
          x.col(eci) = arma::mean(ngb_cols, 1);
          y_est.col(eci) = arma::mean(ngb_cols, 1);
        }
      }
    }
  }
  
  arma::vec rmse_vec(nb);
  
  if (rmse) {
    rmse_vec = rmse_cpp(x, y_est, foc_ind);  
  }

  // compute the number of segments
  arma::uword nseg = cpt.n_elem + 1;
  
  // compute the number of gaps
  arma::uword ngap = ny - ncc;
  
  // assign the year to the gap occurrence
  tpav = tpav % yrs;
  
  // compute the indices at which every segment starts
  arma::uvec seg_beg = arma::join_vert(arma::zeros<arma::uvec>(1), cpt);
  
  // compute length of segments
  arma::vec len_seg(ny, arma::fill::value(arma::datum::nan));

  for (arma::uword i = 0; i < seg_beg.n_elem; i++) {
    if (seg_beg.n_elem == 1) {
      len_seg[seg_beg[i]] = ny;
    }
    else if (i == seg_beg.n_elem - 1) {
      len_seg[seg_beg[i]] = ny - seg_beg[i];
    }
    else {
      len_seg[seg_beg[i]] = seg_beg[i + 1] - seg_beg[i];
    }
  }

  // compute slope values and place them at the beginning of each segment
  arma::mat slo_seg(nb, ny, arma::fill::value(arma::datum::nan));

  for (arma::uword i = 0; i < seg_beg.n_elem; i++) {
    if (len_seg[seg_beg[i]] == 1) {
      slo_seg.col(seg_beg[i]) = arma::zeros(nb);
    }
    else {
      arma::uvec seg_beg_0 = { seg_beg[i] + 1 };
      arma::uvec seg_beg_1 = { seg_beg[i] };
      slo_seg.col(seg_beg[i]) = y_est.submat(foc_ind, seg_beg_0) - y_est.submat(foc_ind, seg_beg_1);
    }
  }
  
  // compute the magnitude of change and identify the type of change
  arma::mat mag(nr, ny, arma::fill::value(arma::datum::nan)), rel_mag(nr, ny, arma::fill::value(arma::datum::nan));
  arma::vec cpt_id(ny, arma::fill::value(arma::datum::nan)), cpt_foc(ny, arma::fill::value(arma::datum::nan));
  arma::rowvec rel_mag_med(ny, arma::fill::value(arma::datum::nan));
  
  double d_max_yr = arma::datum::nan, d_max_mg = arma::datum::nan, 
    d_frs_yr = arma::datum::nan, d_frs_mg = arma::datum::nan, 
    d_lst_yr = arma::datum::nan, d_lst_mg = arma::datum::nan, d_num = arma::datum::nan, 
    g_max_yr = arma::datum::nan, g_max_mg = arma::datum::nan, 
    g_frs_yr = arma::datum::nan, g_frs_mg = arma::datum::nan, 
    g_lst_yr = arma::datum::nan, g_lst_mg = arma::datum::nan, g_num = arma::datum::nan;
  
  if (nseg > 1) {
    arma::uword hnb = floor(nb/2);
    arma::uword hnc = foc_ind[0];
    arma::uvec mat_ind = arma::regspace<arma::uvec>(0, nr - 1);
    
    for (arma::uword i = 0; i < cpt.n_elem; i++) {
      arma::uword cpti = cpt[i];
      arma::umat cptind_cell(nc, nb), cnd_d_mat(nc, nb), cnd_g_mat(nc, nb);
      
      cptind_cell(mat_ind) = cptind.col(i);
      arma::uvec cptind_sum = arma::sum(cptind_cell, 1);
      cpt_foc[cpti] = cptind_sum[hnc];
      
      arma::mat mag1_act(nc, nb), mag1_est(nc, nb);
      mag1_act(mat_ind) = x.col(cpti - 1) - x.col(cpti);
      mag1_est(mat_ind) = y_est.col(cpti - 1) - y_est.col(cpti);
      
      arma::vec cng_dir_act = arma::sign(mag1_act).as_col();
      arma::vec cng_dir_est = arma::sign(mag1_est).as_col();
      
      cnd_d_mat(mat_ind) = cng_dir_act == cng_dir && cng_dir_est == cng_dir;
      cnd_g_mat(mat_ind) = cng_dir_act == -cng_dir && cng_dir_est == -cng_dir;
      
      arma::uvec cnd_d = arma::sum(cnd_d_mat, 1);
      arma::uvec cnd_g = arma::sum(cnd_g_mat, 1);
      
      if (cpti == ny && use_last == 0) {
        cpt_id[cpti] = 9;
      }
      else {  
        if (cnd_d[hnc] >= hnb) {
          cpt_id[cpti] = 101;    // disturbance
          mag.col(cpti) = mag1_est.as_col();
          rel_mag.col(cpti) = (mag1_est.as_col() / y_est.col(cpti - 1)) % arma::vec(nr, arma::fill::value(100));
        }
        else if (cnd_g[hnc] >= hnb) {    
          cpt_id[cpti] = 102;    // greening
          mag.col(cpti) = mag1_est.as_col();
          rel_mag.col(cpti) = (mag1_est.as_col() / y_est.col(cpti - 1)) % arma::vec(nr, arma::fill::value(100));
        }
        else {
          cpt_id[cpti] = 9;
        }
      }
    }
    
    // median change magnitude among bands and cells
    for (arma::uword i = 0; i < cpt.n_elem; i++) {
      arma::vec rel_mag_ci = rel_mag.col(cpt[i]);
      arma::uvec fin_ind = arma::find_finite(rel_mag_ci);
      
      if (fin_ind.n_elem > 0) {
        rel_mag_med[cpt[i]] = arma::median(arma::abs(rel_mag_ci(fin_ind)));
      }
      else {
        continue;
      }
    }

    arma::uvec d_ind = arma::find(cpt_id == 101);
    arma::uvec g_ind = arma::find(cpt_id == 102);
    
    arma::rowvec d_mag_med = rel_mag_med, g_mag_med = rel_mag_med;
    d_mag_med(g_ind).fill(arma::datum::nan);   // remove greening events
    g_mag_med(d_ind).fill(arma::datum::nan);   // remove disturbance events
    
    if (d_ind.n_elem > 0) {
      
      // find the highest-magnitude disturbance and its duration
      arma::uword d_max_ind = arma::index_max(d_mag_med);
      d_max_yr = yrs[d_max_ind];
      d_max_mg = rel_mag_med[d_max_ind];

      // find the first disturbance occurred
      d_frs_yr = yrs[d_ind[0]];
      d_frs_mg = rel_mag_med[d_ind[0]];

      // find the last disturbance occurred
      d_lst_yr = yrs[d_ind[d_ind.n_elem - 1]];
      d_lst_mg = rel_mag_med[d_ind[d_ind.n_elem - 1]];

      // number of disturbances
      d_num = d_ind.n_elem;
    }
    
    if (g_ind.n_elem > 0) {
      
      // find the highest-magnitude greening and its duration
      arma::uword g_max_ind = arma::index_max(g_mag_med);
      g_max_yr = yrs[g_max_ind];
      g_max_mg = rel_mag_med[g_max_ind];

      // find the first greening occurred
      g_frs_yr = yrs[g_ind[0]];
      g_frs_mg = rel_mag_med[g_ind[0]];

      // find the lat greening occurred
      g_lst_yr = yrs[g_ind[g_ind.n_elem - 1]];
      g_lst_mg = rel_mag_med[g_ind[g_ind.n_elem - 1]];

      // number of greening events
      g_num = g_ind.n_elem;
    }
  }

  arma::vec sng_val = {d_max_yr, d_max_mg, d_frs_yr, d_frs_mg, d_lst_yr, d_lst_mg, d_num,
                       g_max_yr, g_max_mg, g_frs_yr, g_frs_mg, g_lst_yr, g_lst_mg, g_num, 
                       static_cast<double>(ngap), static_cast<double>(ntpa)};
  
  arma::vec vec_0 = arma::join_vert(len_seg, cpt_id, cpt_foc, tpav);
  arma::vec vec_val = arma::join_vert(rmse_vec, vec_0, rel_mag_med.t());
  
  arma::vec mat_val = arma::join_vert(y_est.rows(foc_ind).as_col(), slo_seg.as_col(),
                                      mag.rows(foc_ind).as_col(), rel_mag.rows(foc_ind).as_col());
  
  arma::vec res = arma::join_vert(sng_val, vec_val, mat_val);
  return res;
}