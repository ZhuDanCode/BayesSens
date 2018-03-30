// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix kronecker_sp_3_cpp(NumericMatrix B, NumericMatrix A) {
  int step = A.nrow() / B.ncol();
  NumericMatrix C(B.nrow() * step, A.ncol());
  for (int i = 0; i < B.nrow(); i++) {
    for (int j = 0; j < A.ncol(); j++) {
      for (int k = 0; k < B.ncol(); k++) {
        for (int s = 0; s < step; s++) {
          C(step * i + s, j) = C(step * i + s, j) + B(i, k) * A(k * step + s, j);
        }
      }
    }
  }
  return C;
}


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}
