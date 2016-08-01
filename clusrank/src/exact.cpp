#include <Rcpp.h>
using namespace Rcpp;
Rcpp::IntegerVector score_, w, csize;
int meanrks;
//[[Rcpp::export]]

int crksum(int rks, int I, int J, int sumrks, int minrks, int maxrks) {
  int i, j, csize_cum, res;
  Rcpp::IntegerVector score_sub, idx,  score_sum, csize_ord, csize_sub;
  /* Rcpp::Rcout <<rks <<" " <<  I << " " << J  << " " << sumrks << " " << minrks<< std::endl; */
  if (I < 0 | J < 0) return(0);    
  if (I > J) {
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
  Rcpp::Rcout << rks << " " << i << " " << j << " " << sumrks << " " << minrks << " " << maxrks << " " << std::endl;
  if ((rks < minrks) |  (rks > (maxrks))) return(0);
  if (rks == minrks) return(1);
  if (i == 0) {
    return( rks == 0);
  }

  if (j > 0 ) {
    csize_ord = Rcpp::clone(csize);
    std::nth_element(csize_ord.begin(), csize_ord.begin() + i + j, csize_ord.end());
    std::sort(csize_ord.begin(), csize_ord.begin() + i + j);
    csize_sub = csize_ord[seq_len(j) - 1];
    /*   Rcpp::Rcout << "minrks = " << minrks << " rks = " << rks << " i = " << i << " j = " << j <<  std::endl;
        Rcpp::Rcout << "rks -mrks = " << rks - meanrks << " sum(csize_sub) = " << sum(csize_sub) << std::endl; */
    if(((rks - meanrks) < sum(csize_sub))) {
    csize_cum = 0;
    int count = 0;
    for( int ct = 0; csize_cum  < (rks - meanrks); ct = ct + 1) {
      csize_cum += csize_ord[ct];
      count++;
    }
    j = count - 1;
    Rcpp::Rcout << " j = ? " << j << std::endl;
    if ( j < 0 ) return(0);
    if ( j == 0 ) {
      idx = seq_len(i) - 1;
      score_sub = score_[idx];
      int temp;
      temp = sum(score_sub);
      Rcpp::Rcout << "j == 0" << std::endl;
      return(rks == temp);
    }
    idx = seq_len(i + j) - 1;
    score_sum = score_[idx];
    sumrks = sum(score_sum);
    idx = seq_len(i) - 1 + j;
    Rcpp::Rcout <<" count = "<< count <<" idx = " << idx <<" i= " << i << " j = " << j << std::endl; 
        
    score_sub = score_[idx];
    maxrks = sum(score_sub);
    /*  Rcpp::Rcout << rks << " " << i << " " << j << " " << sumrks << " " << minrks << " " << maxrks << " " << std::endl; */
    return(crksum(rks, i, j, sumrks, minrks, maxrks));
    }
  
    
  }
  
  /*   if (rks < j) {
       idx_sum = seq_len(i + rks) - 1;
       score_sum = score_[idx_sum];
       sumrks = sum(score_sum);
       return(crksum(rks, i, rks, sumrks, minrks));
       } 
  */
  if (j == 0) {
    return(rks == sumrks);
  } else {
    sumrks = sumrks - score_[i + j - 1];
    res = crksum(rks - score_[i+j-1], i-1, j, sumrks, minrks - score_[i - 1], maxrks - score_[i + j - 1]) + crksum(rks, i, j-1, sumrks, minrks, maxrks - score_[i + j - 1] + score_[i - 1]);
    return(res);
  }
}

//[[Rcpp::export]]
double pcrksum(int rks, int I, IntegerVector Score, IntegerVector CSize) {
  int sumrks, minrks, maxrks, n, J;
  IntegerVector idx, score_sub;
  int N;
  double  nrksum = 0;
  meanrks = I * (I + 1) / 2;
  score_ = Score;
  csize = CSize;
  n = Score.size();
  J = n - I;
  N = Rf_choose(I + J, I);
  sumrks = sum(Score);
  idx = seq_len(I) - 1;
  score_sub = score_[idx];
  minrks = sum(score_sub); /* Smallest possible rank sum */
  idx = idx + J;
  score_sub = score_[idx];
  maxrks = sum(score_sub); /* Largest possible rank sum */
  for (int k = meanrks; k <= rks; k++) {
    Rcpp::Rcout << "minrks = " << minrks << " rks = " << rks << " I = " << I << " J = " << J <<" k = " << k << " sumrks = " << sumrks <<    std::endl;
    nrksum += crksum(k, I, J, sumrks, minrks, maxrks);
  }
  return(nrksum);
}

//[[Rcpp::export]]
IntegerMatrix cumcrksum(int rks, int I, IntegerVector Score, IntegerVector CSize) {
  IntegerMatrix res(rks + 1, 2);
  int sumrks, minrks, maxrks, n, J;
  IntegerVector idx, score_sub;

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
  
  for (int i = 0; i <= rks; i ++) {
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
       if (zero == 0 & temp <= k) {
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
  for (int i = 0; i < d; i ++) {
    N *= Rf_choose(n[i], xn[i]);
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
  int N, max_s, sum_s, u, j;
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
  for (j = 2; j < max_s + 1; ++j ) {
    compare[0] = j * (j + 1) / 2;
    int i, end = min(compare);
    for (i = end; i >= j; --i){
      w[i] += w[i - j];
      if (i == j & w1[j] == 1) w[i] += 1;
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
