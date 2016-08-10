// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;

template <typename T>
void checkSize(const Map<T> Q1, const Map<T> Q2) {
  if(Q1.size() != Q2.size()) {
    stop("Different vector size");
  }
}

template <typename T>
double ComputeRmse(const Map<T> Q1, const Map<T> Q2) {
  checkSize(Q1, Q2);
  double aux = 0.0;
  int n = Q1.size();
  for (int i = 0; i < Q1.size(); i++) {
    if (R_IsNA(Q1(i)) || R_IsNA(Q2(i))) {
      n--;
    } else {
      aux += (Q1(i) - Q2(i)) * (Q1(i) - Q2(i));
    }
  }
  if (n == 0) {
    return NA_REAL;
  }
  return sqrt(aux / double(n));
}

//' Compute root squared mean error
//' @param Q1 Numeric.
//' @param Q2 Numeric.
//' @export
// [[Rcpp::export]]
double ComputeRmse(SEXP Q1, SEXP Q2) {

  if (TYPEOF(Q1) == REALSXP && TYPEOF(Q2) == REALSXP) {
    return ComputeRmse(as< Map<VectorXd> >(Q1), as< Map<VectorXd> >(Q2));
  } else if (TYPEOF(Q1) == INTSXP && TYPEOF(Q2) == INTSXP) {
    return ComputeRmse(as< Map<VectorXi> >(Q1), as< Map<VectorXi> >(Q2));
  } else {
    stop("Only integer and numeric type are supported, and the vectors must be of same type");
  }
  return NA_REAL;
}



template <typename T>
double ComputeAveragedCrossEntropy(const Map<T> P, const Map<T> Q) {
  checkSize(P, Q);
  double aux = 0.0;
  int n = P.size();
  for (int i = 0; i < P.size(); i++) {
    if (R_IsNA(Q(i)) || R_IsNA(P(i))) {
      n--;
    } else if (Q(i) <= 0) {
      aux += - P(i) * log(1e-6);
    } else {
      aux += -P(i) * log(Q(i));
    }
  }
  if (n == 0) {
    return NA_REAL;
  }
  return aux / double(n);
}

//' Compute average cross entropy
//' @param P Numeric.
//' @param Q Numeric.
//' @export
// [[Rcpp::export]]
double ComputeAveragedCrossEntropy(SEXP P, SEXP Q) {
  if (TYPEOF(P) == REALSXP && TYPEOF(Q) == REALSXP) {
    return ComputeAveragedCrossEntropy(as< Map<VectorXd> >(P), as< Map<VectorXd> >(Q));
  } else if (TYPEOF(P) == INTSXP && TYPEOF(Q) == INTSXP) {
    return ComputeAveragedCrossEntropy(as< Map<VectorXi> >(P), as< Map<VectorXi> >(Q));
  } else {
    stop("Only integer and numeric type are supported, and the vectors must be of same type");
  }
  return NA_REAL;
}
