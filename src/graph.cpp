// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

//' TODO
// [[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeHeatKernelWeight(const Eigen::Map<Eigen::MatrixXd> coord, double sigma) {
        SparseMatrix<double> W(coord.rows(), coord.rows());

        W.reserve(VectorXi::Constant(coord.rows(),0.1 * coord.rows())); // we reserve room for 0.1 % of  non-zeros values per column

        double aux = 0.0;
        for(int i = 0; i < coord.rows(); i++) {
                for(int j = 0; j <=i; j++) { // because W is symetric
                        // aux = exp(- d^2 / sigma^2)
                        aux = std::exp(-(coord.row(i) - coord.row(j)).squaredNorm() / sigma / sigma);
                        if(aux > 1e-10) {
                                W.insert(i,j) = aux;
                        }
                }
        }
        W = SparseMatrix<double>(W.selfadjointView<Lower>());
        W.makeCompressed();
        return W;
}
