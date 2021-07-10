// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
  SEXP ProdMat(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}
