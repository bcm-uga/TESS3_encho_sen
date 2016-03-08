// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                 // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision

//' TODO
// [[Rcpp::export]]
Eigen::MatrixXi ComputeXBin(const Eigen::Map<Eigen::MatrixXi> M, int d) {
        MatrixXi MBin = MatrixXi::Zero(M.rows(), M.cols() * (d+1));
        for (int i = 0; i < M.rows(); i++) {
                for (int l = 0; l < M.cols(); l++) {
                        for (int j = 0; j <= d; j++) {
                                MBin(i, (d+1) * l + j) = (M(i,l) == j);
                        }
                }
        }
        return MBin;
}

//' TODO
// [[Rcpp::export]]
Eigen::MatrixXi ComputeXFromXBin(const Eigen::Map<Eigen::MatrixXi> MBin, int d) {
        MatrixXi M = MatrixXi::Zero(MBin.rows(), MBin.cols() / (d+1));

        // aux matrix
        MatrixXi aux(d+1,1);
        for (int i = 0; i < d+1; i++) {
                aux(i,0) = i;
        }

        for (int i = 0; i < M.rows(); i++) {
                for (int l = 0; l < M.cols(); l++) {
                        M(i, l) = (MBin.block(i,(d+1) * l,1,d+1) * aux)(0,0);
                }
        }
        return M;
}
