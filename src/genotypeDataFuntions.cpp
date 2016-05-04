// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;

//' TODO
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd X2XBin(const Rcpp::NumericMatrix & X, int ploidy) {
        MatrixXd XBin = MatrixXd::Zero(X.rows(), X.cols() * (ploidy + 1));
#ifdef _OPENMP
          #pragma omp parallel for
#endif
        for (int i = 0; i < X.rows(); i++) {
                for (int l = 0; l < X.cols(); l++) {
                        if (NumericMatrix::is_na(X(i,l))) {
                                // missing value
                                XBin.block(i, (ploidy + 1) * l, 1, ploidy + 1) = MatrixXd::Constant(1, ploidy + 1, NA_REAL);
                        } else {
                                for (int j = 0; j <= ploidy; j++) {
                                        XBin(i, (ploidy + 1) * l + j) = (X(i, l) == j);
                                }
                        }
                }
        }
        return XBin;
}

//' TODO
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd XBin2X(const Eigen::Map<Eigen::MatrixXd> XBin, int ploidy) {
        MatrixXd X = MatrixXd::Zero(XBin.rows(), XBin.cols() / (ploidy + 1));

        // aux matrix
        MatrixXd aux(ploidy + 1,1);
        for (int i = 0; i < ploidy + 1; i++) {
                aux(i,0) = i;
        }

        for (int i = 0; i < X.rows(); i++) {
                for (int l = 0; l < X.cols(); l++) {
                        X(i, l) = (XBin.block(i, (ploidy + 1) * l, 1, ploidy + 1) * aux)(0,0);
                }
        }
        return X;
}
