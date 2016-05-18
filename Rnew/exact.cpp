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
int csrk(int srk, int I) {

  if(I < 0) return(0);
  if(srk < 0) return(0);
  if( I == 0)  return(1);            
  if(I > 0) {
    return(csrk(srk - score_[I - 1], I - 1) + csrk(srk, I - 1));
  }
}


//[[Rcpp::export]]
int csrkg(int srk, IntegerVector Score) {
  int N, sum_s, u, c, j;
  IntegerVector compare(2);
  N = Score.size();
  sum_s = sum(Score);
  
  u = sum_s * (sum_s + 1) / 2;
  c = (u / 2);
  if( srk < 0 || srk > u)
    return 0;
  if( w[0] == 1) {
    IntegerVector subw(w.begin(), w.begin() + srk + 1);
    return(sum(subw));
  }
  
  compare[1] = c;
  w = seq_len(sum_s + 1);
  std::fill(w.begin(), w.end(), 0);
  w[Score] = 1;
  w[0] = 1;
  Rcpp::Rcout<<w<<std::endl;
  for( j = 2; j < sum_s + 1; ++j ) {
    compare[0] = j * (j + 1) / 2;
    int i, end = min(compare);
    for( i = end; i > j; --i)
      w[i] += w[i - j];
  }

  IntegerVector subw(w.begin(), w.begin() + srk + 1);  
  return(sum(subw));

}
