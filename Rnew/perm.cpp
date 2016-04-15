// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix perm2 (int N) {
  int n = N;
  int n2 = pow(2, N);
  int i, j;
  NumericMatrix x(n2, n), x2((n2 / 2), n);
  if(n == 2) {
    x( 0, _ ) = NumericVector::create(1, 1);
    x( 1, _) = NumericVector::create(1, -1);
    x( 2, _) = NumericVector::create(-1, 1);
    x( 3, _) = NumericVector::create(-1, -1);
    return x; 
  } else {
    for(i = 0; i < (n2 / 2); i++) {
      x(i, (n - 1)) = 1;
      x((i + n2 / 2), (n - 1)) = -1;
    }
    x2 = perm2(n - 1);
    for(i = 0; i < (n2 / 2); i++) {
      for(j = 0; j < (n - 1); j++) {
        x(i, j) = x2(i, j);
        x((i + n2 / 2), j) = x2(i, j);
      }
    }
    return x;
   }
}



