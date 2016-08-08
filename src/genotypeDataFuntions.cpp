// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;

//' Binary coding of the genotype data matrix.
//' @param X Genotype matrix.
//' @param ploidy Numbet of chromosomes.
//' @param XBin Binary genotype matrix to fill.
//' @export
// [[Rcpp::export]]
void X2XBin(const Rcpp::NumericMatrix & X, int ploidy, Eigen::Map<Eigen::MatrixXd> & XBin) {
        if (XBin.rows() != X.rows() || XBin.cols() != X.cols() * (ploidy + 1)) {
          stop("XBin must be of size nrow(X) * (ncol(X) * (ploidy + 1))");
        }
        #pragma omp parallel for
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
        return;
}

//' Convert the binary genotype matrix into genotype matrix.
//' @param ploidy Numbet of chromosomes.
//' @param XBin Binary genotype matrix.
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
