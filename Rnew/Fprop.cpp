#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]

SEXP Fprop(SEXP x_, SEXP clus_, SEXP nvec_, SEXP m_, SEXP n_) {
  Rcpp::NumericVector x(x_) ;
  Rcpp::NumericVector clus(clus_) ;
  Rcpp::IntegerVector nvec(nvec_);
  const  int m = Rcpp::as<int>(m_), n = Rcpp::as<int>(n_);
  Rcpp::NumericVector Fj(m), Fprop(n);
  for( int i = 0; i < n; i++ ) {
    for( int j = 0; j < m; j++ ) {
      for( int k = 0; k < n; k++) {
        if(clus[k] == (j + 1)) {
          if(x[k] < x[i]) Fj[j] += 1;
          if(x[k] == x[i]) Fj[j] += 0.5;
        }
      }
      Fj[j] = Fj[j] / nvec[j];
      if((j + 1) != clus[i])
        Fprop[i] += Fj[j];     
 
    }
    std::fill(Fj.begin(), Fj.end(), 0);
  }
  return Fprop;
}

//[[Rcpp::export]]

double Fi(double X, int I, NumericVector x, IntegerVector clus, IntegerVector nvec, int N) {
  /* Compute the Fi for X_ in each cluster  */
  /* X_: a value, i_: no. of cluster, x_: all scores */
  /* clus_: cluster id, nvec_: number of objects in each cluster */
  /* n_: total number of obs  */
  
  double Fi = 0;
  for(int i = 0;  i < N;i ++) {
    if(clus[i] == I) {
      if(x[i] < X) {
        Fi += 1;
      }
      if(x[i] == X) {
        Fi += 0.5;
      }
    }
  }
  Fi = Fi / nvec[I - 1];
  return((Fi));
}
      
  
    
//[[Rcpp::export]]

double Ftot(double X,  NumericVector x, IntegerVector clus,
            IntegerVector nvec, int N, int M) {
  /* m_: tot no of cluster */
  double st = 0;
  for( int i = 1; i <= M; i ++) {
    st += Fi(X, i, x, clus, nvec, N);
  }
  return((st));
}

//[[Rcpp::export]]
NumericVector Fivec(NumericVector x,
                    IntegerVector clus, IntegerVector nvec,
                     int N, int M) {
  Rcpp::NumericVector  Fivec(N);
  for( int i = 0; i < N; i ++ ) {
    for( int j = 0; j < M; j ++ ) {
      Fivec[i] += Fi(x[i], (j + 1), x, clus, nvec, N);
    }
  }
  return Fivec;
}

//[[Rcpp::export]]
double Fcom(double X, NumericVector x, IntegerVector clus, IntegerVector nvec,
                   int N, int M) {
  double Fcom = 0;
  for( int i = 0; i < M; i ++ ) {
    Fcom += Fi(X, (i + 1), x, clus, nvec, N) * nvec[i];
  }
  return Fcom / N;
}
