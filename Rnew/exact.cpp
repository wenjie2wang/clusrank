#include <Rcpp.h>
using namespace Rcpp;
Rcpp::IntegerVector score_, w;
int crksum(int rks, int I, int J) {
  int minrks, i, j;
  Rcpp::IntegerVector score_sub, idx = seq_len(I) - 1;
  if( I < 0 | J < 0) return(0);    
  if( rks < 0) {
    if( I >= 0) {
      return(0);
    }
  }
  score_sub = score_[idx];
  minrks = sum(score_sub);
  if(rks < J) return(crksum(rks, I, rks));
  if(rks < minrks) return(0);  
  if( I == 0) return(1);
  if( J == 0) {
    if( sum(score_sub) > rks) return(0);
    return(1);
  }
  if( I > 0 & J >= 0) {
    return(crksum(rks - score_[I+J-1], I-1, J) + crksum(rks, I, J-1));
  }
}

//[[Rcpp::export]]
double pcrksumg(int rks, int I, int J, IntegerVector Score) {
  int sumrks, i;
  int N;
  double p, nrksum;
  score_ = Score;
  N = Rf_choose(I + J, I);
  Rcpp::Rcout<<N<<std::endl;
  sumrks = sum(Score);
  if(rks > sumrks / 2) {
    rks = sumrks - rks - 1;
    i = I;
    I = J;
    J = i;
    nrksum = crksum(rks, I, J);
    return( 1 - nrksum / N );
  } else {
    nrksum = crksum(rks, I, J);    
    return(nrksum / N);
  }
}



//[[Rcpp::export]]

int csrkg(int srk, IntegerVector Score) {
  /* Count the no of combination with a sum rank less than srk */
  /* The sum rank is sum(rank[ rank > 0]) in the programme*/
  /* The input sum rank is sum(sign(rank) * rank) */
  int N, max_s, sum_s,u, c, j;
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

