#include <Rcpp.h>
using namespace Rcpp;
Rcpp::IntegerVector score_, w;
//[[Rcpp::export]]

int crksum(int rks, int I, int J, int sumrks) {
  int minrks, maxrks,  i, j, c;
  Rcpp::IntegerVector score_sub, idx, idx_sum, score_sum;
  if( I < 0 | J < 0) return(0);    

  c = (int) (sumrks / 2); 
  if(rks > c | rks > c) {
    rks = sumrks - rks;
    i = J;
    j = I;
  } else {
    i = I;
    j = J;
  }
  
  idx = seq_len(i) - 1;
  score_sub = score_[idx];
  minrks = sum(score_sub); /* Smallest possible rank sum */
  idx = I + J - idx - 1;
  score_sub = score_[idx];
  maxrks = sum(score_sub);
  
  if(rks < minrks |  (rks > (maxrks))) return(0);

  /* Rcpp::Rcout << " i = " << i << " j = " << j <<"rks = " << rks <<  std::endl; */
  if( i == 0) {
    return( rks == 0);
  }
  
 
  if(rks < j) {
    idx_sum = seq_len(i + rks) - 1;
    score_sum = score_[idx_sum];
    sumrks = sum(score_sum);
    return(crksum(rks, i, rks, sumrks));
  }
  if( j == 0) {

    return( rks == 0);
  } else {
    /*Rcpp::Rcout << rks <<  " score_[i + j - 1] = " << score_[i + j - 1] << std::endl;*/
    sumrks = sumrks - score_[i + j - 1];
    return(crksum(rks - score_[i+j-1], i-1, j, sumrks) + crksum(rks, i, j-1, sumrks));
  }
}

//[[Rcpp::export]]
double pcrksumg(int rks, int I, IntegerVector Score) {
  int sumrks, i, n, J;
  int N;
  double p, nrksum;
  score_ = Score;
  n = Score.size();
  J = n - I;
  N = Rf_choose(I + J, I);
  Rcpp::Rcout<<N<<std::endl;
  sumrks = sum(Score);
  if(rks > sumrks / 2) {
    rks = sumrks - rks - 1;
    i = I;
    I = J;
    J = i;
    nrksum = crksum(rks, I, J, sumrks);
    Rcpp::Rcout << nrksum << std::endl;
    return( 1 - nrksum / N );
  } else {
    nrksum = crksum(rks, I, J, sumrks);
    
    Rcpp::Rcout << nrksum << std::endl;
    return(nrksum / N);
  }
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
  srk = (sum_s - srk) / 2 + srk; 
  u = max_s * (max_s + 1) / 2;
  c = (u / 2);
  if( srk < 0 || srk > u)
    return 0; 

  
  compare[1] = c;
  w = IntegerVector(c + 1);
  w1 = IntegerVector(c + 1);
  w[Score[0]] = 1;
  w[0] = 0;

  w1[Score] = 1;
  Rcpp::Rcout<<w<<std::endl;

  for( j = 2; j < N + 1; ++j ) {
    compare[0] = j * (j + 1) / 2;
    int i, end = min(compare);
    for( i = end; i >= j; --i){
        w[i] += w[i - j];
        if(i == j & w1[j] == 1) w[i] += 1;
        /*  Rcpp::Rcout << "j = " << j << " i = " << i << " w = " << w[i] << std::endl;*/
    }
  }
  w[0] = 1;
  Rcpp::Rcout<<w<<std::endl;
  IntegerVector subw(w.begin(), w.begin() + srk + 1);  
  return(sum(subw));

}

