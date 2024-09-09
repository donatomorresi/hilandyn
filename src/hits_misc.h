#ifndef hits_misc_H
#define hits_misc_H

arma::mat wgt_det_cpp(arma::mat x, arma::umat ind_mat, arma::uword nc, arma::uword nb);

arma::vec match_cpp(arma::uvec a, arma::uvec b);

arma::vec filter_bts_cpp(arma::vec a);

arma::mat orth_matrix_cpp(arma::rowvec d);

arma::vec updating_cpp(arma::umat ee, arma::vec wgt_const, arma::vec wgt_lin, arma::rowvec bts_c, arma::uvec idx);

arma::mat balance_np_cpp(arma::uvec paired, arma::umat ee, arma::uvec idx, arma::uword no_of_current_steps, arma::uword n);

arma::mat balance_p_cpp(arma::umat pr, arma::umat ee_p1, arma::umat ee_p2, arma::uvec idx, arma::uword n);

arma::uvec finding_cp_cpp(Rcpp::List bts_obj, arma::uword n);

Rcpp::List L12_cpp(arma::uword l);

#endif