#include <RcppArmadillo.h>
#include "hits_misc.h"


// [[Rcpp::export]]
arma::vec match_cpp(arma::uvec a, arma::uvec b) {

  arma::vec res(a.n_elem, arma::fill::value(arma::datum::nan));

  for (arma::uword i = 0; i < a.n_elem; ++i) {
    arma::uvec logi_vec = (b == a[i]);
    
    if (any(logi_vec)) {
      arma::uvec ind = arma::find(logi_vec);
      res[i] = arma::as_scalar(ind[0]);
    }
  }
  return res;
}


// [[Rcpp::export]]
arma::vec filter_bts_cpp(arma::vec a) {
  
  double a5a1_a4a2_dif = a[4] * a[0] - a[3] * a[1]; 
  double a2a6_a3a5_dif = a[1] * a[5] - a[2] * a[4]; 
  double a4a3_a1a6_dif = a[3] * a[2] - a[0] * a[5];
  
  double w = - std::sqrt(std::pow(a5a1_a4a2_dif, 2) / (std::pow(a2a6_a3a5_dif, 2) + std::pow(a4a3_a1a6_dif, 2) + std::pow(a5a1_a4a2_dif, 2)));
  double u = w * a2a6_a3a5_dif / a5a1_a4a2_dif;
  double v = w * a4a3_a1a6_dif / a5a1_a4a2_dif;
  
  arma::vec df = {u, v, w};
  
  if (df.has_nan()) {
    
    arma::vec a2 = arma::reverse(a);
    
    a5a1_a4a2_dif = a2[4] * a2[0] - a2[3] * a2[1];
    a2a6_a3a5_dif = a2[1] * a2[5] - a2[2] * a2[4];
    a4a3_a1a6_dif = a2[3] * a2[2] - a2[0] * a2[5];
    
    w = - std::sqrt(std::pow(a5a1_a4a2_dif, 2) / (std::pow(a2a6_a3a5_dif, 2) + std::pow(a4a3_a1a6_dif, 2) + std::pow(a5a1_a4a2_dif, 2)));
    u = w * a2a6_a3a5_dif / a5a1_a4a2_dif;
    v = w * a4a3_a1a6_dif / a5a1_a4a2_dif;
    
    df = {w, v, u};
    
    return df;
  }
  
  return df;
}


// [[Rcpp::export]]
arma::mat orth_matrix_cpp(arma::rowvec d) {
  
  arma::mat33 M;
  M.row(0) = d / std::sqrt(arma::accu(arma::pow(d, 2)));
  
  double u = M(0, 0);
  double v = M(0, 1);
  double w = M(0, 2);
  
  arma::rowvec v1 = {1 - std::pow(u, 2), -u*v, -u*w};
  arma::rowvec v2 = {0, -w, v};
  
  M.row(1) = v1 / std::sqrt(arma::accu(arma::pow(v1, 2)));
  M.row(2) = v2 / std::sqrt(arma::accu(arma::pow(v2, 2)));
  
  return M;
}


// [[Rcpp::export]]
Rcpp::List updating_cpp(arma::umat ee, arma::vec wgt_const, arma::vec wgt_lin, arma::rowvec bts_c, arma::uvec idx) {
  
  arma::mat wc0, wl0, tc0, wc1, wl1, tc1;
  wc0 = wl0 = tc0 = wc1 = wl1 = tc1 = arma::mat(ee.n_rows, 3, arma::fill::none);
  
  for (arma::uword i = 0; i < 3; ++i) {
    
    wc0.col(i) = wgt_const(ee.col(i));
    wl0.col(i) = wgt_lin(ee.col(i));
    tc0.col(i) = bts_c(ee.col(i));
  }
  
  arma::mat m = arma::join_horiz(wc0, wl0);
  arma::mat h(3, m.n_rows, arma::fill::none);
  
  for (arma::uword i = 0; i < m.n_rows; ++i) {
    h.col(i) = filter_bts_cpp(m.row(i).t());
  }
  
  arma::mat M0(9, m.n_rows, arma::fill::none);
  
  for (arma::uword i = 0; i < m.n_rows; ++i) {
    M0.col(i) = orth_matrix_cpp(h.col(i).t()).as_col(); 
  }
  
  //arma::inplace_trans(M0);
  
  arma::umat m_i = { {0,1,2}, {3,4,5}, {6,7,8} };
  
  for(arma::uword i = 0; i < 3; ++i) {
    
    // wc1.col(i) = arma::sum(wc0 % M0.cols(m_i.col(i)), 1);
    // wl1.col(i) = arma::sum(wl0 % M0.cols(m_i.col(i)), 1);
    // tc1.col(i) = arma::sum(tc0 % M0.cols(m_i.col(i)), 1);
    wc1.col(i) = arma::sum(wc0.t() % M0.rows(m_i.col(i)), 0);
    wl1.col(i) = arma::sum(wl0.t() % M0.rows(m_i.col(i)), 0);
    tc1.col(i) = arma::sum(tc0.t() % M0.rows(m_i.col(i)), 0);
  }
  
  arma::uvec eating_up0 = ee.col(0);
  arma::uvec eating_up1 = ee.col(1);
  arma::uvec eaten_up = ee.col(2);
  
  arma::uvec idx1 = idx(arma::find_nonfinite(match_cpp(idx, eaten_up)));
  
  // 1) updating X
  bts_c(eaten_up) = tc1.col(0);
  bts_c(eating_up0) = tc1.col(1);
  bts_c(eating_up1) = tc1.col(2);
  
  // 2) updating weight_const
  wgt_const(eaten_up) = wc1.col(0);
  wgt_const(eating_up0) = wc1.col(1);
  wgt_const(eating_up1) = wc1.col(2);
  
  // 3) updating weight_lin
  wgt_lin(eaten_up) = wl1.col(0);
  wgt_lin(eating_up0) = wl1.col(1);
  wgt_lin(eating_up1) = wl1.col(2);
  
  return Rcpp::List::create(Rcpp::Named("weights_const") = wgt_const,
                            Rcpp::Named("weights_lin") = wgt_lin,
                            Rcpp::Named("bts_coeffs") = bts_c,
                            Rcpp::Named("idx") = idx1,
                            Rcpp::Named("h") = h,
                            Rcpp::Named("tc1") = tc1);
  
}


// [[Rcpp::export]]
arma::mat balance_np_cpp(arma::uvec paired, arma::umat ee, arma::uvec idx, arma::uword no_of_current_steps, arma::uword n) {
  
  arma::umat prd(ee.n_rows, 3);
  
  arma::vec vm = match_cpp(vectorise(ee), paired);
  prd(find_nonfinite(vm)).fill(0);
  prd(find_finite(vm)).fill(1);
  
  arma::mat blnc(3, no_of_current_steps);
  //blnc.fill(arma::datum::nan);
  
  arma::uvec f2 = arma::find(arma::sum(prd.cols(0,1), 1) == 2);
  arma::uvec l2 = arma::find(arma::sum(prd.cols(1,2), 1) == 2);
  arma::uvec f2_l2 = arma::join_vert(f2, l2);
  
  arma::uvec nopair = arma::regspace<arma::uvec>(0, ee.n_rows - 1);
  
  arma::uvec col1 = {0};
  arma::uvec col2 = {1};
  arma::uvec col3 = {2};
  //arma::uvec row12 = {0, 1};
  
  arma::vec ee_l2_col2 = arma::conv_to<arma::vec>::from(ee(l2, col2));
  
  if (f2_l2.n_elem > 0) {
    nopair.shed_rows(f2_l2);
  } 
  
  if (f2.n_elem > 0) {
    
    arma::uvec prtn0 = ee(f2, col3) - ee(f2, col1);
    arma::vec prtn = arma::conv_to<arma::vec>::from(prtn0);
    
    arma::vec prtn3 = prtn + 1;
    arma::vec prtn1 = prtn / prtn3;
    arma::vec prtn2 = 1 / prtn3;
    
    blnc.cols(f2) = arma::join_vert(prtn1, prtn2, prtn3);  // arma::join_horiz(prtn1, prtn2, prtn3).as_col();
    // blnc(row12, f2) = arma::join_horiz(prtn1, prtn2).t();
  }
  
  if (l2.n_elem > 0) {
    
    arma::vec vm = match_cpp(ee(l2, col3), idx);
    arma::vec prtn(vm.n_elem);
    
    for (arma::uword i = 0; i < vm.n_elem; ++i) {
      
      if (std::isnan(vm(i))) {
        
        double n_ee_na = n - ee(i, 2) + 1;
        prtn(i) = n_ee_na;
      }
      else {
        arma::uword vmi = vm(i);
        
        if (idx.in_range(vmi + 1)) {    
          double idx_vmi = idx(vmi + 1);
          prtn(i) = idx_vmi - ee_l2_col2(i);
        }
        else {
          double n_ee_na = n - ee(i, 2) + 1;
          prtn(i) = n_ee_na;
        }
      }
    }
    
    arma::vec prtn3 = prtn + 1;
    arma::vec prtn1 = 1 / prtn3;
    arma::vec prtn2 = prtn / prtn3;
    
    blnc.cols(l2) = arma::join_vert(prtn1, prtn2, prtn3);   //arma::join_horiz(prtn1, prtn2, prtn3).as_col();
    //blnc(row12, l2) = arma::join_horiz(prtn1, prtn2).t();
  }
  
  if (!nopair.is_empty()) {
    
    double onethird = 1;
    onethird /= 3;
    
    blnc.cols(nopair).fill(onethird);
  }
  
  return blnc;
}


// [[Rcpp::export]]
arma::mat balance_p_cpp(arma::umat pr, arma::umat ee_p1, arma::umat ee_p2, arma::uvec idx, arma::uword n) {
  
  arma::mat blnc(3, pr.n_cols);
  //blnc.fill(arma::datum::nan);
  
  arma::vec c1 = arma::conv_to<arma::vec>::from(ee_p1.col(2) - ee_p1.col(0));
  arma::vec c2(ee_p2.n_rows);
  
  arma::vec vm = match_cpp(ee_p2.col(2), idx);
  
  arma::uvec col3 = {2};
  
  for (arma::uword i = 0; i < ee_p2.n_rows; ++i) {
    
    if (std::isnan(vm(i))) {
      
      double n_ee_p1 = n - ee_p1(i, 2);   // + 1
      c2(i) = n_ee_p1;
    }
    else {
      arma::uword vmi = vm(i);
      
      if (idx.in_range(vmi + 1)) {        
        double idx_vmi = idx(vmi + 1);     
        c2(i) = idx_vmi - ee_p1(i, 2);
      }
      else {
        double n_ee_p1 = n - ee_p1(i, 2);   // + 1
        c2(i) = n_ee_p1;
      }
    }
  }
  
  arma::vec c1_c2 = c1 + c2;
  
  blnc.row(0) = c1 / c1_c2;
  blnc.row(1) = c2 / c1_c2;
  blnc.row(2) = c1_c2;
  
  return blnc;
}


// [[Rcpp::export]]
arma::uvec finding_cp_cpp(Rcpp::List bts_obj) {   //Rcpp::NumericVector
  
  const arma::uword n = Rcpp::as<arma::uword>(bts_obj["n"]);
  arma::uvec sameboat = Rcpp::as<arma::uvec>(bts_obj["sameboat"]);
  arma::cube decomp_hist = Rcpp::as<arma::cube>(bts_obj["decomp_hist"]);
  
  arma::umat all_edges(1, 1);
  arma::umat edges = arma::conv_to<arma::umat>::from(arma::vectorise(decomp_hist.row(0)));
  
  if (sameboat.n_elem == 1) {
    all_edges = arma::join_vert(edges, sameboat);
  }
  else {
    all_edges = arma::join_vert(arma::reshape(edges, 3, n - 2), sameboat.t());
  }
  
  arma::umat survived_edges(0, 1);
  
  if (all_edges.n_cols > 1) {
    survived_edges = all_edges.cols(arma::find(arma::abs(arma::vectorise(decomp_hist.tube(2, 0))) > arma::datum::eps));
  }
  
  arma::uvec row_13 = {0, 2};
  arma::uvec row_12 = {0, 1};
  arma::uvec row_123 = {0, 1, 2};
  arma::uvec iv = {0};
  
  arma::umat cp(0, 1);
  
  if ((survived_edges.n_elem > 0) && (survived_edges.n_cols > 1)) {
    arma::uword i = 0;
    
    while (i < survived_edges.n_cols - 1) {
      
      arma::mat part_obj0 = decomp_hist.tube(0, 0, arma::size(1, 2));
      arma::umat part_obj1 = arma::conv_to<arma::umat>::from(part_obj0);
      arma::uvec survived_edges0 = arma::vectorise(survived_edges(1, i, arma::size(2, 1)));
      arma::uvec matched(part_obj1.n_cols);
      
      for (arma::uword j = 0; j < part_obj1.n_cols; j++) {
        
        if (match_cpp(part_obj1.col(j), survived_edges0).is_finite()) {
          matched(j) = 1;
        }
      }
      
      matched = arma::find(matched > 0);
      
      arma::uvec survived_edges_sub1 = arma::diff(arma::vectorise(survived_edges(0, i, arma::size(3, 1))));
      arma::uvec survived_edges_sub2 = arma::diff(arma::vectorise(survived_edges(0, i + 1, arma::size(3, 1))));
      
      arma::uvec iv = {i};
      
      if ((survived_edges(3, i) != 0) && (survived_edges_sub1(0) == 1) && (survived_edges_sub2(0) == 1)) {
        cp = arma::join_vert(cp, survived_edges(row_13, iv));
        i = i + 2;
      }
      else if ((survived_edges(3, i) != 0) && (survived_edges_sub1(1) == 1) && (survived_edges_sub2(0) == 1)) {
        cp = arma::join_vert(cp, survived_edges(row_13, iv + 1));
        i = i + 2;
      }
      else if ((survived_edges(3, i) == 0) && (survived_edges_sub1(0) == 1) && (survived_edges_sub1(1) != 1)) {
        cp = arma::join_vert(cp, survived_edges(row_13, iv));
        i = i + 1;
      }
      else if ((survived_edges(3, i) == 0) && (survived_edges_sub1(0) == 1) && (survived_edges_sub1(1) == 1) && (matched.n_elem > 0)) {
        cp = arma::join_vert(cp, survived_edges(row_12, iv));
        i = i + 1;
      }
      else {
        cp = arma::join_vert(cp, survived_edges(row_123, iv));
        i = i + 1;
      }
    }
    
    cp = arma::unique(cp);
    cp.resize(cp.n_elem + 1);
    cp(cp.n_elem - 1) = n + 1;
  }
  else if ((survived_edges.n_elem > 0)) {
    arma::uvec survived_edges_sub0 = arma::diff(survived_edges(0, 0, arma::size(3, 1)));
    
    if ((survived_edges.n_cols == 1) && (survived_edges.row(3)(0) == 0) && (survived_edges_sub0(0) == 1) && (survived_edges_sub0(1) != 1)) {
      
      cp = arma::join_vert(cp, survived_edges(row_13, iv));
    }
    else if ((survived_edges.n_cols == 1) && (survived_edges.row(3)(0) == 0) && (survived_edges_sub0(0) == 1) && (survived_edges_sub0(1) == 1)) {
      
      cp = arma::join_vert(cp, survived_edges(row_12, iv));
    }
  }
  else {
    arma::umat cp(0, 1);
  }
  
  if ((n == 3) && (cp.n_elem > 0) && (survived_edges.n_cols == 1)) {
    cp.shed_row(0);
    cp.resize(cp.n_elem + 1);
    cp(cp.n_elem - 1) = n;
  }
  else {
    cp = cp(arma::find((cp <= n) && (cp > 1)));
  }
  
  if (cp.n_elem > 0) {
    cp -= 1;
  }
  
  //return Rcpp::NumericVector(cp.begin(), cp.end());
  return cp;
}


// [[Rcpp::export]]
Rcpp::List L12_cpp(arma::uword l) {
  
  arma::uvec x = arma::regspace<arma::uvec>(0, 1, l);
  arma::uword n = x.n_elem;
  
  arma::umat edges(n - 2, 3);
  edges.col(0) = arma::regspace<arma::uvec>(n - 3, 0);            
  edges.col(1) = arma::regspace<arma::uvec>(n - 2, 1);              
  edges.col(2) = arma::regspace<arma::uvec>(n - 1, 2); 
  
  arma::vec weights_const = arma::ones<arma::vec>(n);
  arma::vec weights_lin = arma::regspace<arma::vec>(1, n);
  arma::mat updatedS = arma::diagmat(weights_const);
  
  Rcpp::List L1L2(edges.n_rows);
  
  for (arma::uword i = 0; i < edges.n_rows; i++) {
    
    arma::uvec ee = edges.row(i).as_col();
    arma::vec h = filter_bts_cpp(arma::join_vert(weights_const(ee), weights_lin(ee)));
    arma::mat M = orth_matrix_cpp(h.t());
    
    arma::mat tmp(3, 2);
    tmp.col(0) = weights_const(ee);
    tmp.col(1) = weights_lin(ee);
    arma::mat sm_det = M * tmp;
    
    updatedS.col(ee(0)) = updatedS.cols(edges.row(i)) * M.row(1).t();
    updatedS.col(ee(1)) = updatedS.cols(edges.row(i)) * M.row(2).t();
    updatedS.col(ee(2)) = updatedS.cols(edges.row(i)) * M.row(0).t();
    
    L1L2[i] = updatedS(arma::span(ee(0), n-1), arma::span(ee(0), ee(1)));
    
    arma::uword eating_up0 = ee(0);
    arma::uword eating_up1 = ee(1);
    arma::uword eaten_up = ee(2);
    
    weights_const(eaten_up) = sm_det(0, 0);
    weights_const(eating_up0) = sm_det(1, 0);
    weights_const(eating_up1) = sm_det(2, 0);
    
    weights_lin(eaten_up) = sm_det(0, 1);
    weights_lin(eating_up0) = sm_det(1, 1);
    weights_lin(eating_up1) = sm_det(2, 1);
  }
  
  return L1L2;
}
