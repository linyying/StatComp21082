#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param thin the number of between-sample random numbers
//' @param N the length of chain
//' @param a parameter 
//' @param b parameter
//' @param n the number of samples
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' gibbsC <- gibbs_sampleC(2000,1,4,50,10)
//' qqplot(gibbsC[,1],gibbsR[,1],pch=16,cex=0.7,main="compare x") 
//' abline(a=0,b=1,col = "blue")
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbs_sampleC(int N, double a, double b, int n, int thin) {
  NumericMatrix mat(N, 2);
  double x = floor(n/2), y = 0.5;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y )[0];
      y = rbeta(1, x+a, n-x+b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
