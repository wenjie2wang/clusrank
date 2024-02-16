//
// clusrank: Wilcoxon Rank Tests for Clustered Data
//
// Copyright (C) 2015-2024 Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
// Copyright (C) 2022-2024 Wenjie Wang
//
// This file is part of the R package clusrank.
//
// The R package clusrank is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package clusrank is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//



#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]

NumericVector Fprop(NumericVector x, NumericVector clus,
		    IntegerVector nvec, int m , int n) {

  Rcpp::NumericVector Fj(m), Fprop(n);
  for (int i = 0; i < n; i++ ) {
    for (int j = 0; j < m; j++ ) {
      for (int k = 0; k < n; k++) {
        if (clus[k] == (j + 1)) {
          if (x[k] < x[i]) Fj[j] += 1;
          if (x[k] == x[i]) Fj[j] += 0.5;
        }
      }
      Fj[j] = Fj[j] / nvec[j];
      if ((j + 1) != clus[i])
        Fprop[i] += Fj[j];

    }
    std::fill(Fj.begin(), Fj.end(), 0);
  }
  return Fprop;
}

//[[Rcpp::export]]

double Fi(double X, int I, NumericVector x, IntegerVector clus,
	  IntegerVector nvec, int N) {
  /* Compute the Fi for X_ in each cluster  */
  /* X_: a value, i_: no. of cluster, x_: all scores */
  /* clus_: cluster id, nvec_: number of objects in each cluster */
  /* n_: total number of obs  */

  double Fi = 0;
  for (int i = 0;  i < N;i ++) {
    if (clus[i] == I) {
      if (x[i] < X) {
        Fi += 1;
      }
      if (x[i] == X) {
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
  for (int i = 1; i <= M; i ++) {
    st += Fi(X, i, x, clus, nvec, N);
  }
  return((st));
}

//[[Rcpp::export]]
NumericVector Ftot_vec (NumericVector x, IntegerVector clus,
			IntegerVector nvec, int N, int M) {
  NumericVector res(N);
  for(int i = 0; i < N; i ++ ) {
    res[i] = Ftot(x[i], x, clus, nvec, N, M);
  }
  return(res);
}

//[[Rcpp::export]]
NumericVector Fi_vec(NumericVector x, IntegerVector clus,
		     IntegerVector nvec, int N, int M) {
  Rcpp::NumericVector  Fivec(N);
  for (int i = 0; i < N; i ++ ) {
      Fivec[i] = Fi(x[i], clus[i], x, clus, nvec, N);
  }
  return Fivec;
}

//[[Rcpp::export]]
double Fcom(double X, NumericVector x, IntegerVector clus,
	    IntegerVector nvec, int N, int M) {
  double Fcom = 0;
  for (int i = 0; i < M; i ++ ) {
    Fcom += Fi(X, (i + 1), x, clus, nvec, N) * nvec[i];
  }
  return Fcom / N;
}

//[[Rcpp::export]]
NumericVector Fcom_vec(NumericVector x, IntegerVector clus,
		       IntegerVector nvec, int N, int M) {
  NumericVector result(N);
  for (int i = 0; i < N; i ++) {
    result[i] = Fcom(x[i], x, clus, nvec, N, M);
  }
  return result;
}
