#include <Rcpp.h>
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


// NumericMatrix kronecker_sp_2_cpp(NumericMatrix B, NumericMatrix A) {
//   int nc = B.ncol();
//   int nr = B.nrow();
//   int n = A.nrow() / B.ncol();
//   NumericMatrix C(B.nrow() * n, A.ncol());
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < nr; j++) {
//       for (int k = 0; k < nc; k++) {
//         for (int s = 0; s < ; s++) {
//           C[(1+(i-1)*nr):(i*nr), ] = B(j,k) %*% A[(1+(i-1)*nc):(i*nc), , drop = F];
//         }
//       }
//     }
//   }
//   C
// }
