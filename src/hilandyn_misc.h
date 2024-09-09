#ifndef hilandyn_misc_H
#define hilandyn_misc_H

arma::umat rle_cpp(arma::uvec x);

double mad_cpp(arma::vec x);

arma::vec sd_est_cpp(arma::mat x, arma::vec mad_bias);

arma::vec wgt_sam_cpp(arma::mat x, arma::uword nb, arma::uword nc);

arma::mat focal_values_cpp(arma::vec d, arma::uvec dim, arma::uword w);

arma::vec rmse_cpp(arma::mat x, arma::mat y, arma::uvec foc_ind);

arma::umat cpt_cnd_cpp(arma::mat x, arma::mat y, arma::mat nob, arma::uvec cpt, arma::uword nb, arma::uword nc, arma::uword nob_init_min);

arma::mat proc_cpt_cpp(arma::mat x, arma::umat cptind, arma::umat chk_mat, arma::uvec cpt, arma::vec sd, arma::vec mad_bias, arma::vec wgts, arma::uword nc, arma::uword nb, double th_const);

#endif