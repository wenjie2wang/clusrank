//
// clusrank: Wilcoxon Rank Tests for Clustered Data
// Copyright (C) 2016-2021  Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
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
Rcpp::IntegerVector score_, w;
int meanrks, csize;

int crksum(int rks, int I, int J, int sumrks, int minrks, int maxrks) {
  int i, j, res;
  Rcpp::IntegerVector score_sub, idx,  score_sum;
  /*   Rcpp::Rcout <<rks <<" " <<  I << " " << J  << " " << sumrks << " " << minrks << " " << maxrks << std::endl;  */
  if ((I < 0) | (J < 0)) return(0);
  if (I > J) {
    rks = sumrks - rks;
    i = J;
    j = I;
    if( i != 0) {
      idx = seq_len(i) - 1;
      score_sub = score_[idx];
      minrks = sum(score_sub); /* Smallest possible rank sum*/
      idx = idx + j;
      score_sub = score_[idx];
      maxrks = sum(score_sub);
    }
  } else {
    i = I;
    j = J;
  }
  /*  Rcpp::Rcout << rks << " " << i << " " << j << " " << sumrks << " " << minrks << " " << maxrks << " " << std::endl; */
  if ((rks < minrks) | (rks > maxrks)) return(0);
  if ((rks == minrks) | (rks == maxrks)) return(1);
  if (i == 0) {
    return(rks == 0);
  }


       /*   Rcpp::Rcout << "minrks = " << minrks << " rks = " << rks << " i = " << i << " j = " << j <<  std::endl;
        Rcpp::Rcout << "rks -mrks = " << rks - meanrks << " sum(csize_sub) = " << sum(csize_sub) << std::endl; */
    if(((rks -  i * (i + 1) / 2) < (csize * j))) {
      j = (rks - i * (i + 1) / 2) / csize;
      idx = seq_len(i + j) - 1;
      score_sum = score_[idx];
      sumrks = sum(score_sum);
      idx = seq_len(i) - 1 + j;
      /*   Rcpp::Rcout <<" count = "<< count <<" idx = " << idx <<" i= " << i << " j = " << j << std::endl; */
      score_sub = score_[idx];
      maxrks = sum(score_sub);
      return(crksum(rks, i, j, sumrks, minrks, maxrks));
    }

    sumrks = sumrks - score_[i + j - 1];
    // Rcpp::Rcout << "maxrks = " << maxrks << std::endl;
    res = crksum(rks - score_[i+j-1], i-1, j, sumrks, minrks - score_[i - 1], maxrks - score_[i + j - 1]) + crksum(rks, i, j-1, sumrks, minrks, maxrks - score_[i + j - 1] + score_[j - 1]);
    return(res);

}

double pcrksum(int rks, int I, IntegerVector Score, int Csize) {
  int sumrks, minrks, maxrks, n, J;
  IntegerVector idx, score_sub;
  double  nrksum = 0;
  meanrks = I * (I + 1) / 2;
  score_ = Score;
  csize = Csize;
  n = Score.size();
  J = n - I;
  sumrks = sum(Score);
  idx = seq_len(I) - 1;
  score_sub = score_[idx];
  minrks = sum(score_sub); /* Smallest possible rank sum */
  idx = idx + J;
  score_sub = score_[idx];
  maxrks = sum(score_sub); /* Largest possible rank sum */
  for (int k = meanrks; k <= rks; k++) {
    /*  Rcpp::Rcout << "minrks = " << minrks << " rks = " << rks << " I = " << I << " J = " << J <<" k = " << k << " sumrks = " << sumrks <<    std::endl; */
    nrksum += crksum(k, I, J, sumrks, minrks, maxrks);
  }
  return(nrksum);
}

//[[Rcpp::export]]
IntegerMatrix cumcrksum(int rks, int I, IntegerVector Score, int Csize) {
  IntegerMatrix res(rks + 1, 2);
  int sumrks, minrks, maxrks, n, J;
  IntegerVector idx, score_sub;

  csize = Csize;
  meanrks = I * (I + 1) / 2;
  score_ = Score;
  n = Score.size();
  J = n - I;
  sumrks = sum(Score);

  idx = seq_len(I) - 1;
  score_sub = score_[idx];
  minrks = sum(score_sub); /* Smallest possible rank sum */
  idx = idx + J;
  score_sub = score_[idx];
  maxrks = sum(score_sub); /* Largest possible rank sum */

  for (int i = meanrks; i <= rks; i++) {
    res(i, 1) = crksum(i, I, J, sumrks, minrks, maxrks);
    res(i, 0) = i;
  }
  return(res);
}



//[[Rcpp::export]]
int crksum_str(int k, IntegerMatrix x, IntegerMatrix xc, IntegerVector max) {
  /* x: matrix of non zero rank sum, xc: matrix of count of non zero rank sum */
  int d = x.ncol(), temp = 0,  ct = 0, ctp = 1;
  IntegerVector slots(d), y(d);
  int id = 0, zero = 0, counter = 0;
  while(true) {
    counter += 1;
       for (int i = 0; i < d; i ++ ) {
         if(xc(slots[i], i) < 0) {
           zero = 1;
           break;
         }
         temp += x(slots[i], i);
       }
       if ((zero == 0) & (temp <= k)) {
         for (int i = 0; i < d; i ++ ) {
           ctp *= xc(slots[i], i);
         }
         ct += ctp;
       }
       temp = 0;
       ctp = 1;
       slots[0]++;
       while (slots[id] == max[id]) {
         if (id == d - 1) {
           return(ct);
         }
         slots[id++] = 0;
         slots[id]++;
       }
       id = 0;
  }
}


//[[Rcpp::export]]
double pcrksum_str(int k, IntegerMatrix x, IntegerMatrix xc, IntegerVector xn, IntegerVector n, IntegerVector max) {
  int d = x.ncol(), N = 1;
  double p = 0;
  for (int i = 0; i < d; i++) {
    N *= Rf_choose(n[i], xn[i]);
  }
  p = crksum_str(k, x, xc, max);
  p = p / N;
  return(p);
}







int csrkg(int srk, IntegerVector Score) {
  /* Count the no of combination with a sum rank less than srk */
  /* The sum rank is sum(rank[ rank > 0]) in the programme*/
  /* The input sum rank is sum(sign(rank) * rank) */
  int N, max_s, u, j, k;
  IntegerVector compare(2), w1;
  N = Score.size();
  max_s = max(Score);
  u = max_s * (max_s + 1) / 2;
  if((srk < 0) | (srk > u))
    return 0;
  compare[1] = u;
  w = IntegerVector(u + 1);
  w1 = IntegerVector(u + 1);
  w[0] = 0;

  for (int i = 0; i < N; ++ i) {
    k = Score[i];
    w1[k] = w1[k] + 1;
  }


  for (j = 2; j < max_s + 1; ++j ) {
    compare[0] = j * (j + 1) / 2;
    int i, end = min(compare);
    for (i = end; i >= j; --i){

      if (w1[j] == 1) {
	w[i] += w[i - j] ;
      } else if (w1[j] > 1) {
	w[i] += w[i - j] * Rf_choose(w1[j], 1);
      	for ( int k = 2; k <= w1[j]; ++k ) {
	  if (i - k * j > 0) {
	    w[i] += w[i - k * j] * Rf_choose(w1[j], k);
	  } else if (i == k * j) {
	    w[i] += Rf_choose(w1[j], k);
	  }
      	}
      }
      if ((w1[i] != 0) & (i == j)) {
	w[i] += w1[i];
      }

    }
  }
  w[0] = 1;
  IntegerVector subw(w.begin(), w.begin() + srk + 1);
  return(sum(subw));
}

//[[Rcpp::export]]

double psrkg(int srk, IntegerVector Score) {
  int n0 = Score.size();
  double n = (double)n0;
  double N = pow(2.0, n), p = 0;
  p = csrkg(srk, Score);
  p = p / N;
  return(p);
}
