#include <RcppArmadillo.h>
#include "hits_misc.h"


// [[Rcpp::export]]
Rcpp::List hd_bts_dcmp_cpp(arma::mat bts_coeffs, double p, arma::vec w) {
  
  const arma::uword d = bts_coeffs.n_rows;
  const arma::uword n = bts_coeffs.n_cols;
  arma::uword noe = n - 2;
  
  arma::vec weights_const = arma::ones<arma::vec>(n);
  arma::vec weights_lin = arma::regspace<arma::vec>(1, n);
  arma::uvec idx = arma::regspace<arma::uvec>(0, n-1);              // starts from 0
  arma::uvec paired;                                    
  
  arma::umat edges(noe, 4);                                         // filled with zeros by default 
  edges.col(0) = arma::regspace<arma::uvec>(0, n - 3);              // starts from 0
  edges.col(1) = arma::regspace<arma::uvec>(1, n - 2);              // starts from 1
  edges.col(2) = arma::regspace<arma::uvec>(2, n - 1);              // starts from 2
  
  arma::cube decomp_hist(4 * d, 3, n - 2, arma::fill::none);
  
  arma::uword steps_left = n - 2;
  arma::uword current_step = 0;
  arma::uvec sameboat;
  
  arma::uvec ind1 = {0};
  arma::uvec ind2 = {1};
  arma::uvec ind3 = {2};
  arma::uvec ind4 = {3};
  arma::uvec ind12 = {0,1};
  arma::uvec ind258 = {1,4,7};
  arma::uvec ind369 = {2,5,8};
  
  while (edges.n_rows > 0) {
    
    double steps_left_d = steps_left;
    arma::uword max_current_steps = std::ceil(p * steps_left_d);
    arma::uvec removable_nodes(idx.max() + 1, arma::fill::ones);
    
    arma::mat Dmat(d, edges.n_rows);
    
    if (any(edges.col(3) > 0)) {
      
      arma::uvec prv = arma::find(edges.col(3) != 0);
      arma::umat pr = arma::reshape(prv, 2, prv.n_elem / 2);
      
      arma::uvec pr_row1 = pr.row(0).t();
      arma::uvec pr_row2 = pr.row(1).t();
      arma::rowvec weights_const_pr = weights_const(edges(pr_row2, ind3)).t();
      arma::rowvec weights_lin_pr = weights_lin(edges(pr_row2, ind3)).t();

      for (arma::uword j = 0; j < d; ++j) {
        
        arma::uvec row_j = {j};
        
        // start computation of detail coefficients
        arma::mat sub_wc, sub_wl, sub_tc, detcoef;
        
        sub_wc = sub_wl = sub_tc = arma::mat(pr_row1.n_elem, 3, arma::fill::none);
        detcoef = arma::mat(3, pr_row1.n_elem, arma::fill::none);
        
        arma::umat ind = edges.rows(pr_row1);
        
        for (arma::uword i = 0; i < 3; ++i) {
          
          sub_wc.col(i) = weights_const(ind.col(i));
          sub_wl.col(i) = weights_lin(ind.col(i));
          sub_tc.col(i) = bts_coeffs(row_j, ind.col(i)).t();
        }
        
        arma::mat sub_m = arma::join_horiz(sub_wc, sub_wl).t();
        
        for (arma::uword i = 0; i < pr_row1.n_elem; ++i) {
          detcoef.col(i) = filter_bts_cpp(sub_m.col(i));
        }
        
        arma::rowvec details = arma::sum(detcoef % sub_tc.t(), 0);
        // end computation of detail coefficients

        arma::mat M0(9, detcoef.n_cols, arma::fill::none);
        
        for (arma::uword k = 0; k < detcoef.n_cols; ++k) {
          M0.col(k) = orth_matrix_cpp(detcoef.col(k).t()).as_col();
        }
        
        arma::uvec rowj(1);
        rowj[0] = j;

        arma::mat cd1_wc = sub_wc.t();
        arma::mat cd1_wl = sub_wl.t();
        arma::mat cd1_tc = sub_tc.t();
        
        arma::mat upd_wc = arma::join_vert(sum(cd1_wc % M0.rows(ind258), 0),
                                           sum(cd1_wc % M0.rows(ind369), 0),
                                           weights_const_pr);

        arma::mat upd_wl = arma::join_vert(sum(cd1_wl % M0.rows(ind258), 0),
                                           sum(cd1_wl % M0.rows(ind369), 0),
                                           weights_lin_pr);

        arma::mat upd_bts_coeffs = arma::join_vert(sum(cd1_tc % M0.rows(ind258), 0),
                                                   sum(cd1_tc % M0.rows(ind369), 0),
                                                   bts_coeffs(rowj, edges(pr_row2, ind3)));

        arma::mat p2m = arma::join_vert(upd_wc, upd_wl);
        arma::mat p2d(3, pr.n_cols, arma::fill::none);
        
        for (arma::uword q = 0; q < pr.n_cols; ++q) {
          p2d.col(q) = filter_bts_cpp(p2m.col(q));
        }
        
        arma::vec p2_details = sum(p2d % upd_bts_coeffs, 0).t();
        
        arma::vec p_detail0 = arma::max(arma::join_horiz(arma::abs(details.t()), arma::abs(p2_details)), 1);
        
        arma::vec p_detail = arma::repelem(p_detail0, 2, 1);
        
        if (pr.n_elem != edges.n_rows) {
          
          arma::uvec edgerow_cd3 = arma::regspace<arma::uvec>(0, edges.n_rows - 1);
          edgerow_cd3.shed_rows(arma::vectorise(pr));

          // start computation of detail coefficients
          arma::mat sub_wc, sub_wl, sub_tc, detcoef;
          
          sub_wc = sub_wl = sub_tc = arma::mat(edgerow_cd3.n_elem, 3, arma::fill::none);
          detcoef = arma::mat(3, edgerow_cd3.n_elem, arma::fill::none);
          
          arma::umat ind = edges.rows(edgerow_cd3);
          
          for (arma::uword i = 0; i < 3; ++i) {
            sub_wc.col(i) = weights_const(ind.col(i));
            sub_wl.col(i) = weights_lin(ind.col(i));
            sub_tc.col(i) = bts_coeffs(row_j, ind.col(i)).t();
          }
          
          arma::mat sub_m = arma::join_horiz(sub_wc, sub_wl).t();
          
          for (arma::uword i = 0; i < edgerow_cd3.n_elem; ++i) {
            detcoef.col(i) = filter_bts_cpp(sub_m.col(i));
          }
          
          arma::rowvec details = arma::sum(detcoef % sub_tc.t(), 0);
          // end computation of detail coefficients
          
          Dmat.row(j) = arma::join_vert(p_detail, details.t()).t();
        }
        else {
          Dmat.row(j) = p_detail.t();
        }
      }
    }
    else {
      
      arma::uvec edgerow = arma::regspace<arma::uvec>(0, edges.n_rows - 1);
      arma::mat sub_wc = arma::join_horiz(weights_const(edges(edgerow, ind1)), weights_const(edges(edgerow, ind2)), weights_const(edges(edgerow, ind3)));
      arma::mat sub_wl = arma::join_horiz(weights_lin(edges(edgerow, ind1)), weights_lin(edges(edgerow, ind2)), weights_lin(edges(edgerow, ind3)));
      arma::mat sub_m = arma::join_horiz(sub_wc, sub_wl).t();
      
      arma::mat detcoef(3, sub_m.n_cols, arma::fill::none); 

      for (arma::uword q = 0; q < detcoef.n_cols; ++q) {
        detcoef.col(q) = filter_bts_cpp(sub_m.col(q));
      }

      Dmat = bts_coeffs.cols(edges.col(0)) % arma::repmat(detcoef.row(0), bts_coeffs.n_rows, 1) + 
        bts_coeffs.cols(edges.col(1)) % arma::repmat(detcoef.row(1), bts_coeffs.n_rows, 1) + 
        bts_coeffs.cols(edges.col(2)) % arma::repmat(detcoef.row(2), bts_coeffs.n_rows, 1);
    }

    arma::mat Dmat2 = arma::abs(Dmat) % arma::repmat(w, 1, Dmat.n_cols);
  
    arma::rowvec colmaxD = arma::max(Dmat2, 0) + arma::mean(Dmat2, 0); 
  
    arma::uvec ord_det = arma::stable_sort_index(colmaxD);

    arma::umat cand = arma::join_horiz(ord_det, edges(ord_det, ind4));
    
    arma::uvec eitr(1);
    arma::uword tei = 0;
    
    if (cand(0, 1) > 0) {
      
      removable_nodes(edges(ord_det.subvec(0,1), ind1)).zeros();
      removable_nodes(edges(ord_det.subvec(0,1), ind2)).zeros();
      removable_nodes(edges(ord_det.subvec(0,1), ind3)).zeros();
      tei += 1;
      eitr = arma::resize(eitr, eitr.n_elem + 1, 1);
      eitr.tail(1) = tei;
    }
    else {
      removable_nodes(edges(ord_det(0), 0)) = removable_nodes(edges(ord_det(0), 1)) = removable_nodes(edges(ord_det(0), 2)) = 0;
    }

    while ((eitr.n_elem < max_current_steps) & (tei < noe)) {
      
      tei += 1;
      
      if (cand(tei, 1) > 0) {
        
        arma::uvec edges_ind1 = edges(ord_det.subvec(tei, tei + 1), ind1);
        arma::uvec edges_ind2 = edges(ord_det.subvec(tei, tei + 1), ind2);
        arma::uvec edges_ind3 = edges(ord_det.subvec(tei, tei + 1), ind3);
        
        if ((arma::sum(removable_nodes(edges_ind1)) == 2) && 
            (arma::sum(removable_nodes(edges_ind2)) == 2) && 
            (arma::sum(removable_nodes(edges_ind3)) == 2)) {
          
          removable_nodes(edges_ind1).zeros();
          removable_nodes(edges_ind2).zeros();
          removable_nodes(edges_ind3).zeros();
          
          eitr = arma::resize(eitr, eitr.n_elem + 2, 1);
          eitr.tail(2) = {tei, tei + 1};
          tei += 1;
        }
      }
      else {
        
        arma::uvec edges_ind1(1), edges_ind2(1), edges_ind3(1);
        
        edges_ind1 = edges(ord_det(tei), 0);
        edges_ind2 = edges(ord_det(tei), 1);
        edges_ind3 = edges(ord_det(tei), 2);
        
        if (!removable_nodes(edges_ind1).is_empty() && !removable_nodes(edges_ind2).is_empty() && !removable_nodes(edges_ind3).is_empty()) {
          
          eitr = arma::resize(eitr, eitr.n_elem + 1, 1);
          eitr.tail(1) = tei;
          removable_nodes(arma::join_vert(edges_ind1, edges_ind2, edges_ind3)).zeros();
          
        }
      }
    }
    
    arma::uvec details_min_ind = ord_det(eitr);
    
    arma::uword no_of_current_steps = eitr.n_elem;
    
    arma::umat ee = arma::reshape(edges.rows(details_min_ind), no_of_current_steps, 4);
    
    arma::uvec idx0(idx.n_elem);
    idx0 = idx;
    
    arma::uvec ee_col4_ind = arma::find(ee.col(3) > 0);
    
    if (ee_col4_ind.is_empty()) {
      
      sameboat = arma::join_vert(sameboat, ee.col(3));
      ee.shed_col(3);
      
      Rcpp::List udt(6);
      
      arma::uword slice0 = current_step;
      arma::uword slice1 = current_step + no_of_current_steps - 1;

      for (arma::uword j = 0; j < d; ++j) {
        
        udt = updating_cpp(ee, weights_const, weights_lin, bts_coeffs.row(j), idx0);
        
        bts_coeffs.row(j) = Rcpp::as<arma::rowvec>(udt["bts_coeffs"]);
        
        decomp_hist(arma::span(j * 4), arma::span::all, arma::span(slice0, slice1)) = arma::conv_to<arma::mat>::from(ee.t() + 1);
        decomp_hist(arma::span(j * 4 + 1), arma::span::all, arma::span(slice0, slice1)) = Rcpp::as<arma::mat>(udt["h"]);
        decomp_hist(arma::span(j * 4 + 2), arma::span::all, arma::span(slice0, slice1)) = Rcpp::as<arma::mat>(udt["tc1"]).t();
        decomp_hist(arma::span(j * 4 + 3), arma::span::all, arma::span(slice0, slice1)) = balance_np_cpp(paired, ee, idx0, no_of_current_steps, n);
      }
      
      weights_const = Rcpp::as<arma::vec>(udt["weights_const"]);
      weights_lin = Rcpp::as<arma::vec>(udt["weights_lin"]);
      idx = Rcpp::as<arma::uvec>(udt["idx"]);
    }
    else {
      
      sameboat = arma::join_vert(sameboat, arma::vectorise(ee(arma::find(ee.col(3) != 0), ind4)));
      
      arma::uvec pr0 = arma::find(ee.col(3) != 0);
      arma::umat pr = arma::reshape(pr0, 2, pr0.n_elem / 2);

      ee(pr.row(1), ind12) = ee(pr.row(0), ind12);
      
      arma::umat ee_p1 = ee.rows(pr.row(0));
      arma::umat ee_p2 = ee.rows(pr.row(1));
      ee_p1.shed_col(3); ee_p2.shed_col(3);
      
      Rcpp::List udt(6);

      arma::uword slice0 = current_step;
      arma::uword slice1 = current_step + ee_p1.n_rows - 1;
      arma::uword slice2 = current_step + ee_p1.n_rows;
      arma::uword slice3 = current_step + pr.n_elem - 1;
      arma::uword slice4 = current_step + pr.n_elem;
      arma::uword slice5 = current_step + no_of_current_steps - 1;
      
      for (arma::uword j = 0; j < d; ++j) {
        
        udt = updating_cpp(ee_p1, weights_const, weights_lin, bts_coeffs.row(j), idx0);

        decomp_hist(arma::span(j * 4), arma::span::all, arma::span(slice0, slice1)) = arma::conv_to<arma::mat>::from(ee_p1.t() + 1);
        decomp_hist(arma::span(j * 4 + 1), arma::span::all, arma::span(slice0, slice1)) = Rcpp::as<arma::mat>(udt["h"]);
        decomp_hist(arma::span(j * 4 + 2), arma::span::all, arma::span(slice0, slice1)) = Rcpp::as<arma::mat>(udt["tc1"]).t();
        decomp_hist(arma::span(j * 4 + 3), arma::span::all, arma::span(slice0, slice1)) = balance_p_cpp(pr, ee_p1, ee_p2, idx0, n);

        udt = updating_cpp(ee_p2, Rcpp::as<arma::vec>(udt["weights_const"]), Rcpp::as<arma::vec>(udt["weights_lin"]), Rcpp::as<arma::mat>(udt["bts_coeffs"]), Rcpp::as<arma::uvec>(udt["idx"]));
        idx = Rcpp::as<arma::uvec>(udt["idx"]);
        
        decomp_hist(arma::span(j * 4), arma::span::all, arma::span(slice2, slice3)) = arma::conv_to<arma::mat>::from(ee_p2.t() + 1);
        decomp_hist(arma::span(j * 4 + 1), arma::span::all, arma::span(slice2, slice3)) = Rcpp::as<arma::mat>(udt["h"]);
        decomp_hist(arma::span(j * 4 + 2), arma::span::all, arma::span(slice2, slice3)) = Rcpp::as<arma::mat>(udt["tc1"]).t();
        decomp_hist(arma::span(j * 4 + 3), arma::span::all, arma::span(slice2, slice3)) = balance_p_cpp(pr, ee_p1, ee_p2, idx, n);
        
        if (pr.n_elem != ee.n_rows) {
          
          sameboat = arma::join_vert(sameboat, arma::vectorise(ee(arma::find(ee.col(3) == 0), ind4)));
          arma::umat ee_np = ee; ee_np.shed_rows(arma::vectorise(pr));
          
          udt = updating_cpp(ee_np, Rcpp::as<arma::vec>(udt["weights_const"]), Rcpp::as<arma::vec>(udt["weights_lin"]), Rcpp::as<arma::mat>(udt["bts_coeffs"]), idx);
          bts_coeffs.row(j) = Rcpp::as<arma::mat>(udt["bts_coeffs"]);
          arma::uword no_of_current_steps_pr = no_of_current_steps - pr.n_elem;
          
          decomp_hist(arma::span(j * 4), arma::span::all, arma::span(slice4, slice5)) = arma::conv_to<arma::mat>::from(ee_np.t() + 1);
          decomp_hist(arma::span(j * 4 + 1), arma::span::all, arma::span(slice4, slice5)) = Rcpp::as<arma::mat>(udt["h"]);
          decomp_hist(arma::span(j * 4 + 2), arma::span::all, arma::span(slice4, slice5)) = Rcpp::as<arma::mat>(udt["tc1"]).t();
          decomp_hist(arma::span(j * 4 + 3), arma::span::all, arma::span(slice4, slice5)) = balance_np_cpp(paired, ee_np, idx, no_of_current_steps_pr, n);
        }
        else {
          bts_coeffs.row(j) = Rcpp::as<arma::mat>(udt["bts_coeffs"]);
        }
      }
      weights_const = Rcpp::as<arma::vec>(udt["weights_const"]);
      weights_lin = Rcpp::as<arma::vec>(udt["weights_lin"]);
      idx = Rcpp::as<arma::uvec>(udt["idx"]);

      ee.shed_col(3);
    }
    
    ////// STEP 5: Updating other variables //////
    paired = arma::unique(arma::join_vert(paired, arma::vectorise(ee.cols(0,1))));
    
    arma::vec paired_ee = match_cpp(paired, vectorise(ee.col(2)));
    
    if (!arma::find_finite(paired_ee).is_empty()) {
      paired = arma::sort(paired(arma::find_nonfinite(paired_ee)));
    }
    
    arma::umat edges2(idx.n_elem - 2, 3, arma::fill::none);
    
    if (edges2.n_rows == 0) {
      break;
    }
    
    edges = edges2;
    edges.col(0) = idx(arma::regspace<arma::uvec>(0, idx.n_elem - 3));
    edges.col(1) = idx(arma::regspace<arma::uvec>(1, idx.n_elem - 2));
    edges.col(2) = idx(arma::regspace<arma::uvec>(2, idx.n_elem - 1));
    
    arma::vec edges_paired = match_cpp(arma::vectorise(edges), paired);
    edges_paired(arma::find_finite(edges_paired)).ones();
    edges_paired(arma::find_nonfinite(edges_paired)).zeros();
    
    arma::mat matchpair = arma::reshape(edges_paired, edges_paired.n_elem / 3, 3).t();
    arma::uvec rs = arma::find(arma::sum(matchpair, 0) == 3);

    if (rs.n_elem > 0) {
      
      arma::umat edges_rs = edges; 
      edges_rs.shed_rows(rs);
      
      arma::mat matchpair_rs = matchpair; 
      matchpair_rs.shed_cols(rs);

      edges = arma::join_vert(edges.rows(rs), edges_rs);
      matchpair = arma::join_horiz(matchpair.cols(rs), matchpair_rs);

      arma::uvec edges_seq = arma::repelem(arma::regspace<arma::uvec>(1, std::floor(rs.n_elem / 2)), 2, 1); 
      arma::uvec edges_zero(edges.n_rows - rs.n_elem);
      edges = arma::join_horiz(edges, arma::join_vert(edges_seq, edges_zero));
    }
    else {
      arma::uvec edges_zero(edges.n_rows);
      edges = arma::join_horiz(edges, edges_zero);
    }

    arma::uvec removed = arma::join_vert(arma::find(arma::sum(matchpair, 0) == 1), 
                                         arma::find((matchpair.row(0) == 1) && (matchpair.row(1) == 0) && (matchpair.row(2) == 1)));

    if (removed.n_elem > 0) {
      edges.shed_rows(removed);
    }
    
    noe = edges.n_rows;
    steps_left = steps_left - no_of_current_steps;
    current_step = current_step + no_of_current_steps;

    if (noe == 1) {
      edges.col(3) = 0;
    }
  }

  return Rcpp::List::create(Rcpp::Named("n") = n,
                            Rcpp::Named("sameboat") = sameboat.t(),
                            Rcpp::Named("decomp_hist") = decomp_hist,
                            Rcpp::Named("bts_coeffs") = bts_coeffs);
}


// [[Rcpp::export]]
Rcpp::List hd_bts_dns_cpp(Rcpp::List bts_obj, double lambda, double bal, arma::vec w, arma::uvec foc_ind) {

  const arma::uword n = Rcpp::as<arma::uword>(bts_obj["n"]);
  const arma::uvec sameboat = Rcpp::as<arma::uvec>(bts_obj["sameboat"]);
  arma::cube decomp_hist = Rcpp::as<arma::cube>(bts_obj["decomp_hist"]);

  const arma::uword d = decomp_hist.n_rows / 4;
  const arma::uword s = decomp_hist.n_slices;
  arma::uvec slice_ind = arma::regspace<arma::uvec>(1, d) * 4 - 2;
  arma::mat detail_all(d, s);
  arma::vec details(s);
  
  if (d > 1) {

    for (arma::uword i = 0; i < s; i++) {
      detail_all.col(i) = decomp_hist.slice(i)(slice_ind);
    }
    
    double w_sum = arma::sum(w);

    for (arma::uword i = 0; i < s; i++) {

      arma::vec det_i = arma::abs(detail_all.col(i));

      double det_all = arma::sum(det_i % w) / w_sum;
      double det_foc = arma::mean(det_i(foc_ind));

      if (det_foc > det_all) {
        details[i] = det_foc;
      }
      else {
        details[i] = det_all;
      }
    }
  }
  else {
    detail_all = arma::vectorise(decomp_hist.tube(2, 0)).t();
    details = arma::abs(detail_all.as_col());
  }

  arma::uvec protect(n);

  for (arma::uword i = 0; i < (n - 2); i++) {
    
    arma::uword ind_c1 = decomp_hist(0, 0, i) - 1;
    arma::uword ind_c2 = decomp_hist(0, 1, i) - 1;
    arma::uword ind_c3 = decomp_hist(0, 2, i) - 1;
    

    if ((protect(ind_c1) == 0) && (protect(ind_c2) == 0) && (protect(ind_c3) == 0)) {

      arma::uword logi_val = (details(i) > lambda) && (decomp_hist(3, 0, i) > bal) && (decomp_hist(3, 1, i) > bal);
      arma::uvec logi_vec(slice_ind.n_elem, arma::fill::value(logi_val));

      decomp_hist.slice(i)(slice_ind) = decomp_hist.slice(i)(slice_ind) % logi_vec;
    }
    
    if (std::abs(decomp_hist(2, 0, i)) > 0) {     

      protect(ind_c1) = 1;
      protect(ind_c2) = 1;
    }
  }

  arma::uvec paired0 = arma::find(sameboat != 0);
  arma::umat paired = arma::reshape(paired0, 2, paired0.n_elem / 2);

  if (paired.n_elem > 0)  {
    for (arma::uword i = 0; i < paired.n_cols; i++) {

      arma::vec overzero = match_cpp(paired.col(i), arma::find(arma::abs(decomp_hist.tube(2, 0)) > 0));
      arma::vec zero = match_cpp(paired.col(i), arma::find(arma::abs(decomp_hist.tube(2, 0)) == 0));

      arma::uvec iv = { i };
      arma::uvec sl = paired(arma::find_nonfinite(overzero), iv);

      arma::uvec overzero1 = arma::find_finite(overzero);
      arma::uvec zero1 = arma::find_finite(zero);

      if ((overzero1.n_elem == 1) & (zero1.n_elem == 1)) {
        decomp_hist.slice(arma::as_scalar(sl))(slice_ind) = detail_all.col(arma::as_scalar(sl));
      }
    }
  }

  bts_obj["decomp_hist"] = decomp_hist;
  return bts_obj;
}


// [[Rcpp::export]]
Rcpp::List hd_bts_inv_cpp(Rcpp::List bts_obj) {
  
  const arma::uword n = Rcpp::as<arma::uword>(bts_obj["n"]);
  arma::cube decomp_hist = Rcpp::as<arma::cube>(bts_obj["decomp_hist"]);
  arma::mat bts_coeffs = Rcpp::as<arma::mat>(bts_obj["bts_coeffs"]);
  
  const arma::uword d = decomp_hist.n_rows / 4;
  const arma::uvec slice_ind_1 = arma::regspace<arma::uvec>(1, d) * 4 - 2;
  const arma::uvec add_row(slice_ind_1.n_elem, arma::fill::value(decomp_hist.n_rows));
  const arma::uvec slice_ind_2 = slice_ind_1 + add_row;
  const arma::uvec slice_ind_3 = slice_ind_2 + add_row;
  
  arma::mat33 inv_mat;
  
  for (arma::uword i = (n - 2); i-- > 0; ) {
    
    inv_mat = orth_matrix_cpp(decomp_hist(1, 0, i, arma::size(1, 3, 1))).t();

    arma::uvec ind = arma::conv_to<arma::uvec>::from(arma::vectorise(decomp_hist(0, 0, i, arma::size(1, 3, 1)))) - 1;

    decomp_hist.slice(i)(slice_ind_2) = bts_coeffs.col(ind(0));
    decomp_hist.slice(i)(slice_ind_3) = bts_coeffs.col(ind(1));
    
    arma::mat tmp = arma::reshape(decomp_hist.slice(i)(arma::join_vert(slice_ind_1, slice_ind_2, slice_ind_3)), slice_ind_1.n_elem, 3);
    arma::mat rcstr_tmp(1, 1);
    
    if (d == 1) {
      rcstr_tmp = inv_mat * tmp.as_col();
    }
    else {
      rcstr_tmp = inv_mat * tmp.t();
    }

    bts_coeffs.col(ind(0)) = rcstr_tmp.row(0).t();
    bts_coeffs.col(ind(1)) = rcstr_tmp.row(1).t();
    bts_coeffs.col(ind(2)) = rcstr_tmp.row(2).t();
  }
  
  bts_obj["decomp_hist"] = decomp_hist;
  bts_obj["bts_coeffs"] = bts_coeffs;
  return bts_obj;
}


// [[Rcpp::export]]
Rcpp::List hd_bts_pp1_cpp(Rcpp::List bts_obj, double lambda) {

  const arma::uword n = Rcpp::as<arma::uword>(bts_obj["n"]);
  arma::mat pp1fit = Rcpp::as<arma::mat>(bts_obj["bts_coeffs"]);

  arma::rowvec wc = arma::ones<arma::rowvec>(n);
  arma::rowvec wl = arma::regspace<arma::rowvec>(1, n);

  arma::uvec inicp = finding_cp_cpp(bts_obj);
  
  arma::uvec ind_13 = {0, 1, 2};
  
  if (inicp.n_elem > 0) {
    arma::uvec chp = arma::join_vert(arma::ones<arma::uvec>(1),
                                     inicp + arma::ones<arma::uvec>(inicp.n_elem),
                                     arma::uvec(1, arma::fill::value(n + 1)));

    arma::umat pqr = arma::join_horiz(chp(arma::regspace<arma::uvec>(0, chp.n_elem - 3)),
                                      chp(arma::regspace<arma::uvec>(1, chp.n_elem - 2)) - 1,
                                      chp(arma::regspace<arma::uvec>(2, chp.n_elem - 1)) - 1);
    
    pqr -= arma::ones<arma::umat>(arma::size(pqr));

    arma::imat d_pqr = arma::diff(arma::conv_to<arma::imat>::from(pqr), 1, 1);

    arma::mat details, detail_1, detail_2;
    details = detail_1 = detail_2 = arma::mat(pp1fit.n_rows, d_pqr.n_rows);
    
    while (chp.n_elem > 2) {

      arma::uvec c1 = arma::find((d_pqr.col(0) == 0) && (d_pqr.col(1) == 1));
      
      if (c1.n_elem > 0) {

        for (arma::uword i = 0; i < c1.n_elem; i++) {
          arma::uword p = arma::as_scalar(pqr(c1[i], 0));
          arma::uword q = arma::as_scalar(pqr(c1[i], 1));
          arma::uword r = arma::as_scalar(pqr(c1[i], 2));
          details.col(c1[i]) = arma::abs((pp1fit.cols(p, q) - pp1fit.cols(q + 1, r)) / std::sqrt(2));
        }
      }

      arma::umat cnd1 = (d_pqr.col(0) == 0) && (d_pqr.col(1) == 2).as_col();
      arma::umat cnd2 = (d_pqr.col(0) == 1) && (d_pqr.col(1) == 1).as_col();
      arma::uvec c23 = arma::find(cnd1 || cnd2); 
      
      if (c23.n_elem > 0) {

        for (arma::uword i = 0; i < c23.n_elem; i++) {
          arma::uword p = arma::as_scalar(pqr(c23[i], 0));
          arma::uword r = arma::as_scalar(pqr(c23[i], 2));
          arma::vec h = filter_bts_cpp(arma::trans(arma::join_horiz(wc.subvec(p, r), wl.subvec(p, r))));
          details.col(c23[i]) = arma::abs(pp1fit.cols(p, r) * h);
        }
      }

      arma::uvec c4 = arma::find((d_pqr.col(0) == 1) && (d_pqr.col(1) == 2));
      
      if (c4.n_elem > 0) {

        for (arma::uword i = 0; i < c4.n_elem; i++) {
          arma::uword p1 = arma::as_scalar(pqr(c4[i], 0));
          arma::uword q1 = arma::as_scalar(pqr(c4[i], 1));
          arma::uword r1 = q1 + 1;
          arma::vec h1 = filter_bts_cpp(arma::trans(arma::join_horiz(wc.subvec(p1, r1), wl.subvec(p1, r1))));
          detail_1.col(c4[i]) = arma::abs(pp1fit.cols(p1, r1) * h1);
          
          arma::mat M = orth_matrix_cpp(h1.t());
          
          arma::uword r2 = arma::as_scalar(pqr(c4[i], 2, arma::size(1, 1)));
          
          arma::mat M1 = M * wc.subvec(p1, r1).t();
          arma::mat M2 = M * wl.subvec(p1, r1).t();

          arma::vec h2 = filter_bts_cpp(arma::join_vert(M1(1, 0, arma::size(2, 1)), arma::vec(1, arma::fill::value(wc[r2])), 
                                                        M2(1, 0, arma::size(2, 1)), arma::vec(1, arma::fill::value(wl[r2]))));
          
          arma::mat M_pp1 = arma::trans(M * pp1fit.cols(p1, r1).t());
          
          detail_2.col(c4[i]) = arma::abs(arma::join_horiz(M_pp1.cols(1, 2), pp1fit.col(r2)) * h2);
          
          details.col(c4[i]) = arma::max(arma::join_horiz(arma::abs(detail_1.col(c4[i])), arma::abs(detail_2.col(c4[i]))), 1);
          }
        }
      
        arma::uvec c5 = arma::find(d_pqr.col(0) == 0 && d_pqr.col(1) > 2); 
          
          if (c5.n_elem > 0) {

            for (arma::uword i = 0; i < c5.n_elem; i++) {
              arma::uword p = arma::as_scalar(pqr(c5[i], 0));
              arma::uword q = arma::as_scalar(pqr(c5[i], 1));
              arma::uword r = arma::as_scalar(pqr(c5[i], 2));
              
              arma::mat l12 = L12_cpp(r - q)[r - q - 3];
              arma::rowvec wcL = wc.subvec(q + 1, r) * l12;
              arma::rowvec wlL = wl.subvec(q + 1, r) * l12;
              arma::mat xL = pp1fit.cols(q + 1, r) * l12;

              arma::vec h = filter_bts_cpp(arma::trans(arma::join_horiz(wc.subvec(p, q), wcL, wl.subvec(p, q), wlL)));
              details.col(c5[i]) = arma::abs(arma::join_horiz(pp1fit.cols(p, q), xL) * h);
              }
            }
          
          arma::uvec c6 = arma::find(d_pqr.col(0) > 1 && d_pqr.col(1) == 1); 
          
          if (c6.n_elem > 0) {

            for (arma::uword i = 0; i < c6.n_elem; i++) {
              arma::uword p = arma::as_scalar(pqr(c6[i], 0));
              arma::uword q = arma::as_scalar(pqr(c6[i], 1));
              arma::uword r = arma::as_scalar(pqr(c6[i], 2));
              
              arma::mat l12 = L12_cpp(q - p + 1)[q - p - 2];
              arma::rowvec wcL = wc.subvec(p, q) * l12;
              arma::rowvec wlL = wl.subvec(p, q) * l12;
              arma::mat xL = pp1fit.cols(p, q) * l12;
              
              arma::vec h = filter_bts_cpp(arma::trans(arma::join_horiz(wcL, wc.subvec(q + 1, r), wlL, wl.subvec(q + 1, r))));
              details.col(c6[i]) = arma::abs(arma::join_horiz(xL, pp1fit.cols(q + 1, r)) * h);
            }
          }
          
          arma::uvec c7 = arma::find(d_pqr.col(0) == 1 && d_pqr.col(1) > 2);
          
          if (c7.n_elem > 0) {

            for (arma::uword i = 0; i < c7.n_elem; i++) {
              arma::uword p = arma::as_scalar(pqr(c7[i], 0));
              arma::uword q = arma::as_scalar(pqr(c7[i], 1));
              arma::uword r = arma::as_scalar(pqr(c7[i], 2));
              
              arma::mat l12 = L12_cpp(r - q)[r - q - 3];
              arma::rowvec wcL = wc.subvec(q + 1, r) * l12;
              arma::rowvec wlL = wl.subvec(q + 1, r) * l12;
              arma::mat xL = pp1fit.cols(q + 1, r) * l12;
              
              arma::rowvec new_wc = arma::join_horiz(wc.subvec(p, q), wcL);
              arma::rowvec new_wl = arma::join_horiz(wl.subvec(p, q), wlL);
              arma::mat new_x = arma::join_horiz(pp1fit.cols(p, q), xL);
              
              arma::vec h1 = filter_bts_cpp(arma::trans(arma::join_horiz(new_wc.subvec(0, 2), new_wl.subvec(0, 2))));
              detail_1.col(c7[i]) = arma::abs(new_x.cols(0, 2) * h1);
              
              arma::mat M = orth_matrix_cpp(h1.t());
              
              arma::vec M_new_wc = M * new_wc(ind_13);
              arma::vec M_new_wl = M * new_wl(ind_13);
              arma::vec new_wc_4 = arma::vec(1, arma::fill::value(new_wc[3]));
              arma::vec new_wl_4 = arma::vec(1, arma::fill::value(new_wl[3]));
              
              arma::vec h2 = filter_bts_cpp(arma::join_vert(M_new_wc(1, 0, arma::size(2, 1)), new_wc_4,
                                                            M_new_wl(1, 0, arma::size(2, 1)), new_wl_4));
              
              arma::mat M_new_x = arma::trans(M * new_x.cols(0, 2).t());
              detail_2.col(c7[i]) = arma::abs(arma::join_horiz(M_new_x.cols(1, 2), new_x.col(3)) * h2);
              details.col(c7[i]) = arma::max(arma::join_horiz(arma::abs(detail_1.col(c7[i])), arma::abs(detail_2.col(c7[i]))), 1);
            }
          }

          arma::uvec c8 = arma::find(d_pqr.col(0) > 1 && d_pqr.col(1) == 2);
          
          if (c8.n_elem > 0) {

            for (arma::uword i = 0; i < c8.n_elem; i++) {
              arma::uword p = arma::as_scalar(pqr(c8[i], 0));
              arma::uword q = arma::as_scalar(pqr(c8[i], 1));
              arma::uword r = arma::as_scalar(pqr(c8[i], 2));
              
              arma::mat l12 = L12_cpp(q - p + 1)[q - p - 2];
              arma::rowvec wcL = wc.subvec(p, q) * l12;
              arma::rowvec wlL = wl.subvec(p, q) * l12;
              arma::mat xL = pp1fit.cols(p, q) * l12;

              arma::rowvec new_wc = arma::join_horiz(wcL, wc.subvec(q + 1, r));
              arma::rowvec new_wl = arma::join_horiz(wlL, wl.subvec(q + 1, r));
              arma::mat new_x = arma::join_horiz(xL, pp1fit.cols(q + 1, r));
              
              arma::vec h1 = filter_bts_cpp(arma::trans(arma::join_horiz(new_wc.subvec(0, 2), new_wl.subvec(0, 2))));
              detail_1.col(c8[i]) = arma::abs(new_x.cols(0, 2) * h1);
              
              arma::mat M = orth_matrix_cpp(h1.t());
              
              arma::vec M_new_wc = M * new_wc(ind_13);
              arma::vec M_new_wl = M * new_wl(ind_13);
              arma::vec new_wc_4 = arma::vec(1, arma::fill::value(new_wc[3]));
              arma::vec new_wl_4 = arma::vec(1, arma::fill::value(new_wl[3]));
              
              arma::vec h2 = filter_bts_cpp(arma::join_vert(M_new_wc(1, 0, arma::size(2, 1)), new_wc_4,
                                                            M_new_wl(1, 0, arma::size(2, 1)), new_wl_4));
              
              arma::mat M_new_x = arma::trans(M * new_x.cols(0, 2).t());
              detail_2.col(c8[i]) = arma::abs(arma::join_horiz(M_new_x.cols(1, 2), new_x.col(3)) * h2);
              details.col(c8[i]) = arma::max(arma::join_horiz(arma::abs(detail_1.col(c8[i])), arma::abs(detail_2.col(c8[i]))), 1);
            }
          }
          
          arma::uvec all_idx = arma::regspace<arma::uvec>(0, d_pqr.n_rows - 1);
          arma::uvec c1_c8 = arma::join_vert(arma::join_vert(c1, c23, c4, c5), 
                                             arma::join_vert(c6, c7, c8));
          arma::uvec c9 = all_idx(arma::find_nonfinite(match_cpp(all_idx, arma::unique(c1_c8))));
          
          if (c9.n_elem > 0) {

            for (arma::uword i = 0; i < c9.n_elem; i++) {

              arma::uword p = arma::as_scalar(pqr(c9[i], 0));
              arma::uword q = arma::as_scalar(pqr(c9[i], 1));
              arma::uword r = arma::as_scalar(pqr(c9[i], 2));
              
              arma::mat l12_1 = L12_cpp(q - p + 1)[q - p - 2];
              
              arma::mat wcL1 = wc.subvec(p, q) * l12_1;
              arma::mat wlL1 = wl.subvec(p, q) * l12_1;
              arma::mat xL1 = pp1fit.cols(p, q) * l12_1;
              
              arma::mat l12_2 = L12_cpp(r - q)[r - q - 3];
              arma::mat wcL2 = wc.subvec(q + 1, r) * l12_2;
              arma::mat wlL2 = wl.subvec(q + 1, r) * l12_2;
              arma::mat xL2 = pp1fit.cols(q + 1, r) * l12_2;
              
              arma::rowvec new_wc = arma::join_horiz(wcL1, wcL2);
              arma::rowvec new_wl = arma::join_horiz(wlL1, wlL2);
              arma::mat new_x = arma::join_horiz(xL1, xL2);
              
              arma::vec h1 = filter_bts_cpp(arma::trans(arma::join_horiz(new_wc.subvec(0, 2), new_wl.subvec(0, 2))));
              detail_1.col(c9[i]) = arma::abs(new_x.cols(0, 2) * h1);
              
              arma::mat M = orth_matrix_cpp(h1.t());
              arma::vec M_new_wc = M * new_wc(ind_13);
              arma::vec M_new_wl = M * new_wl(ind_13);
              arma::vec new_wc_4 = arma::vec(1, arma::fill::value(new_wc[3]));
              arma::vec new_wl_4 = arma::vec(1, arma::fill::value(new_wl[3]));
              
              arma::vec h2 = filter_bts_cpp(arma::join_vert(M_new_wc(1, 0, arma::size(2, 1)), new_wc_4,
                                                            M_new_wl(1, 0, arma::size(2, 1)), new_wl_4));
              
              arma::mat M_new_x = arma::trans(M * new_x.cols(0, 2).t());
              detail_2.col(c9[i]) = arma::abs(arma::join_horiz(M_new_x.cols(1, 2), new_x.col(3)) * h2);
              details.col(c9[i]) = arma::max(arma::join_horiz(arma::abs(detail_1.col(c9[i])), arma::abs(detail_2.col(c9[i]))), 1);
            }
          }
          
          arma::rowvec maxdet = arma::max(details, 0);
          arma::uword maxdetmin_i = arma::index_min(maxdet);

          if (maxdet[maxdetmin_i] < lambda) {
            
            chp.shed_row(maxdetmin_i + 1);
            
            for (arma::uword k = 0; k < chp.n_elem - 1; k++) {
              
              arma::uvec domain = arma::regspace<arma::uvec>(chp[k] - 1, chp[k + 1] - 2);

              if (domain.n_elem == 1) {
                continue;
              }
              else {
                for (arma::uword q = 0; q < pp1fit.n_rows; q++) {
                  
                  arma::uvec q_vec = {q};
                  arma::vec p1 = arma::polyfit(arma::conv_to<arma::vec>::from(domain), pp1fit(q_vec, domain), 1);
                  pp1fit(q_vec, domain) = arma::trans(arma::polyval(p1, arma::conv_to<arma::vec>::from(domain)));
                }
              }
            }
          }
          else {
            break;
          }
          
          // update
          arma::uvec pi = arma::regspace<arma::uvec>(0, chp.n_elem - 3);
          arma::uvec qi = arma::regspace<arma::uvec>(1, chp.n_elem - 2);
          arma::uvec ri = arma::regspace<arma::uvec>(2, chp.n_elem - 1);
          
          arma::uvec chp_p, chp_q, chp_r;
          
          if (pi.n_elem > 0) {
            chp_p = arma::uvec(pi.n_elem);
            
            for (arma::uword i = 0; i < pi.n_elem; i++) {
              if (chp.in_range(pi[i])) {
                chp_p[i] = chp[pi[i]];
              }
            }
          }

          if (qi.n_elem > 0) {
            chp_q = arma::uvec(qi.n_elem);
            
            for (arma::uword i = 0; i < qi.n_elem; i++) {
              if (chp.in_range(qi[i])) {
                chp_q[i] = chp[qi[i]] - 1;
              }
            }
          }

          if (ri.n_elem > 0) {
            chp_r = arma::uvec(ri.n_elem);
            for (arma::uword i = 0; i < ri.n_elem; i++) {
              if (chp.in_range(ri[i])) {
                chp_r[i] = chp[ri[i]] - 1;
              }
            }
          }

          pqr = arma::join_vert(chp_p, chp_q, chp_r);
          
          arma::uvec nelem = {chp_p.n_elem, chp_q.n_elem, chp_r.n_elem};
          arma::uword nr = arma::max(nelem);
          arma::uword nc = pqr.n_elem / nr;
          pqr = arma::reshape(pqr, nr, nc);
          pqr -= arma::ones<arma::umat>(arma::size(pqr));

          d_pqr = arma::diff(arma::conv_to<arma::imat>::from(pqr), 1, 1);
          details = detail_1 = detail_2 = arma::mat(pp1fit.n_rows, d_pqr.n_rows);
    }
    
    // final change points
    chp.shed_row(0);
    chp.shed_row(chp.n_elem - 1);
    bts_obj["chp"] = chp;
    bts_obj["details"] = details;
    bts_obj["bts_coeffs"] = pp1fit;
    
  }
  else {
    arma::mat details = arma::mat(pp1fit.n_rows, 2);
    bts_obj["chp"] = inicp;
    bts_obj["details"] = details;
    bts_obj["bts_coeffs"] = pp1fit;
  }
  
  return bts_obj;
}


// [[Rcpp::export]]
Rcpp::List hd_bts_cpt_cpp(arma::mat x, arma::vec sd, arma::uword nb, double th_const, arma::vec weights, arma::uvec foc_ind) {

  const double p = 0.01;
  const double bal = 0;
  const arma::uword n = x.n_cols;

  x /= arma::repmat(sd, 1, n);
  double lambda = th_const * std::sqrt(2 * std::log(nb * n));
  
  if (weights.has_nan()) {
    weights = arma::ones(weights.n_elem);
  }

  Rcpp::List dcmp = hd_bts_dcmp_cpp(x, p, weights);
  Rcpp::List dns = hd_bts_dns_cpp(dcmp, arma::as_scalar(lambda), bal, weights, foc_ind);
  Rcpp::List inv = hd_bts_inv_cpp(dns);
  Rcpp::List pp1 = hd_bts_pp1_cpp(inv, arma::as_scalar(lambda));

  arma::umat cptind(arma::size(Rcpp::as<arma::mat>(pp1["details"])));
  cptind(arma::find(Rcpp::as<arma::mat>(pp1["details"]) > arma::as_scalar(lambda))).ones();

  return Rcpp::List::create(Rcpp::Named("est") = Rcpp::as<arma::mat>(pp1["bts_coeffs"]) % arma::repmat(sd, 1, n),
                            Rcpp::Named("cpt") = Rcpp::as<arma::uvec>(pp1["chp"]),
                            Rcpp::Named("cptind") = cptind);
}