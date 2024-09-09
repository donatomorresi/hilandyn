#ifndef hits_H
#define hits_H

Rcpp::List hd_bts_dcmp_cpp(arma::mat bts_coeffs, double p, arma::vec w);

Rcpp::List hd_bts_dns_cpp(Rcpp::List bts_obj, arma::uword n, double lambda, double bal, arma::vec w, arma::uvec foc_ind);

Rcpp::List hd_bts_inv_cpp(Rcpp::List bts_obj, arma::uword n);

Rcpp::List hd_bts_pp1_cpp(Rcpp::List bts_obj, arma::uword n, double lambda);

Rcpp::List hd_bts_cpt_cpp(arma::mat x, arma::vec sd, arma::uword nb, double th_const, arma::vec weights, arma::uvec foc_ind);

#endif