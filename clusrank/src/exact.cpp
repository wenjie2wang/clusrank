#include <Rcpp.h>
using namespace Rcpp;
Rcpp::IntegerVector score_, w;
//[[Rcpp::export]]

int crksum(int rks, int I, int J, int sumrks, int minrks) {
  int i, j, c;
  Rcpp::IntegerVector score_sub, idx, idx_sum, score_sum;
  if( I < 0 | J < 0) return(0);    
  c = (int) (sumrks / 2);
  if(rks > c) {
    rks = sumrks - rks;
    i = J;
    j = I;
    if( i != 0) {
      idx = seq_len(i) - 1;
      score_sub = score_[idx];
      minrks = sum(score_sub); /* Smallest possible rank sum*/
    } else {
      minrks = 0;
    }
  } else {
    i = I;
    j = J;
  }
  
  if(rks < minrks |  (rks > (sumrks))) return(0);
  if( i == 0) {
    return( rks == 0);
  }
  
  
  if(rks < j) {
    idx_sum = seq_len(i + rks) - 1;
    score_sum = score_[idx_sum];
    sumrks = sum(score_sum);
    return(crksum(rks, i, rks, sumrks, minrks));
  }
  if( j == 0) {

    return( rks == 0);
  } else {
    /*Rcpp::Rcout << rks <<  " score_[i + j - 1] = " << score_[i + j - 1] << std::endl;*/
    sumrks = sumrks - score_[i + j - 1];
    
    return(crksum(rks - score_[i+j-1], i-1, j, sumrks, minrks - score_[i - 1]) +
           crksum(rks, i, j-1, sumrks, minrks));
  }
}

//[[Rcpp::export]]
double pcrksum(int rks, int I, IntegerVector Score) {
  int sumrks, minrks, i, n, J;
  IntegerVector idx, score_sub;
  int N;
  double p, nrksum = 0;
  score_ = Score;
  n = Score.size();
  J = n - I;
  N = Rf_choose(I + J, I);
  sumrks = sum(Score);
  
  idx = seq_len(I) - 1;
  score_sub = score_[idx];
  minrks = sum(score_sub); /* Smallest possible rank sum*/
  for( int k = 1; k <= rks; k++) {  
    nrksum += crksum(k, I, J, sumrks, minrks);
  }
  return(nrksum);
}

//[[Rcpp::export]]
IntegerMatrix cumcrksum(int rks, int I, IntegerVector Score) {
  IntegerMatrix res(rks + 1, 2);
  int sumrks, minrks, i, n, J;
  IntegerVector idx, score_sub;

  score_ = Score;
  n = Score.size();
  J = n - I;
  sumrks = sum(Score);
  
  idx = seq_len(I) - 1;
  score_sub = score_[idx];
  minrks = sum(score_sub); /* Smallest possible rank sum*/
  for( int i = 0; i <= rks; i ++) {
    res(i, 1) = crksum(i, I, J, sumrks, minrks);
    res(i, 0) = i;
  }
  return(res);
}



//[[Rcpp::export]]
int crksum_str(int k, IntegerMatrix x, IntegerMatrix xc, IntegerVector max) {
  /* x: matrix of non zero rank sum, xc: matrix of count of non zero rank sum */
  int d = x.ncol(), m = x.nrow(), temp = 0, tempc = 0, ct = 0, ctp = 1;
  IntegerVector slots(d), y(d);
  int id = 0, zero = 0, counter = 0;
  while(true) {
    counter += 1;
       for( int i = 0; i < d; i ++ ) {
         if(xc(slots[i], i) < 0) {
           zero = 1;
           break;
         }
         temp += x(slots[i], i);
       }
       if(zero == 0 & temp < k) {
         for( int i = 0; i < d; i ++ ) {
           ctp *= xc(slots[i], i);
         }
         ct += ctp;
       }
       temp = 0;
       ctp = 1;
       slots[0]++;
       while(slots[id] == max[id]) {
         if( id == d - 1) {
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
  for( int i = 0; i < d; i ++) {
    N *= Rf_choose(n[i], xn[i]);
    Rcpp::Rcout << N << std::endl;

  }
  p = crksum_str(k, x, xc, max);
  p = p / N;
  return(p);
}





//[[Rcpp::export]]

int csrkg(int srk, IntegerVector Score) {
  /* Count the no of combination with a sum rank less than srk */
  /* The sum rank is sum(rank[ rank > 0]) in the programme*/
  /* The input sum rank is sum(sign(rank) * rank) */
  int N, max_s, sum_s,u, c, j
    ;
  IntegerVector compare(2), w1;
  N = Score.size();
  max_s = max(Score);
  sum_s = sum(Score);
  u = max_s * (max_s + 1) / 2;
  if( srk < 0 || srk > u)
    return 0; 
  compare[1] = u;
  w = IntegerVector(u + 1);
  w1 = IntegerVector(u + 1);
  w[Score[0]] = 1;
  w[0] = 0;

  w1[Score] = 1;
  Rcpp::Rcout<<w<<std::endl;

  for( j = 2; j < max_s + 1; ++j ) {
    Rcpp::Rcout << "j = " << j << std::endl;
    compare[0] = j * (j + 1) / 2;
    int i, end = min(compare);
    Rcpp::Rcout << "end = " << end << std::endl;
    for( i = end; i >= j; --i){
      w[i] += w[i - j];
      if(i == j & w1[j] == 1) w[i] += 1;
      /*  Rcpp::Rcout << "j = " << j << " i = " << i << " w = " << w[i] << std::endl;*/
    }
  }
  w[0] = 1;
  Rcpp::Rcout<<w << " srk " <<srk <<std::endl;
  IntegerVector subw(w.begin(), w.begin() + srk + 1);  
  return(sum(subw));
}


//[[Rcpp::export]]

double psrkg(int srk, IntegerVector Score) {
  int n = Score.size();
  double N = pow(n, 2), p = 0;
  p = csrkg(srk, Score);
  p = p / N;
  return(p);
}
