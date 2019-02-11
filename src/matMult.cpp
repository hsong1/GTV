#include <RcppEigen.h>
using namespace Rcpp;

//'@export 
// [[Rcpp::export]]
SEXP matMult(const Eigen::Map<Eigen::MatrixXd> A, 
             const Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

//'@export 
// [[Rcpp::export]]
SEXP matMult_sp(const Eigen::Map<Eigen::SparseMatrix<double> > A, 
             const Eigen::Map<Eigen::SparseMatrix<double> > B){
  Eigen::SparseMatrix<double> C = A * B;
  return Rcpp::wrap(C);
}



/*** R
A <- matrix(rnorm(10000), 100, 100)
B <- matrix(rnorm(10000), 100, 100)
rbenchmark::benchmark(A%*%B,matMult(A,B))

library(Matrix)

Asp = matrix(0, 100,100)
Asp[sample(1:10000,100)] = rnorm(100)
Bsp[sample(1:10000,100)] = rnorm(100)
Asp <- Matrix(Asp, sparse = T)
Bsp <- Matrix(Bsp, sparse = T)

rbenchmark::benchmark(A%*%B,matMult(A,B),Asp%*%Bsp, matMult_sp(Asp,Bsp))

*/
