#ifndef GET_COOCURRENCE_MATRIX_H
#define GET_COOCURRENCE_MATRIX_H
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
IntegerMatrix rcpp_get_coocurrence_matrix(const IntegerMatrix x,
                                          const arma::imat directions,
                                          const int n_cores = 1);

#endif // GET_COOCURRENCE_MATRIX_H
