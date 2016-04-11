// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

using namespace Eigen;

//' TODO
// [[Rcpp::export]]
Eigen::MatrixXd ComputeFst(const Eigen::Map<Eigen::MatrixXd> Q, const Eigen::Map<Eigen::MatrixXd> G, int D) {

        // Fst = 1 - sigma2_s / sigma2_T
        // sigma2_s = Sum_k( q_k * f_k * ( 1 - f_k ) )
        // sigma2_T = ( Sum_k( q_k * f_k ) * ( 1 - Sum_k( q_k * f_k ) ) )

        int K = Q.cols();
        int L = G.rows() / D;

        MatrixXd aux = MatrixXd::Constant(1, D-1, 0.0);
        for (int i = 0; i < aux.rows(); i++) {
          aux(0, i) = double(i+1) / double(D-1);
        }

        // Compute allele frequency
        MatrixXd All(G.rows() / D, K);
        for (int l = 0; l < All.rows(); l++) {
                All.block(l, 0, 1, K) = aux * G.block(D * l + 1, 0, D-1, K);

        }

        //Compute q
        MatrixXd q(1, K);
        for (int k = 0; k < K; k++) {
                q(0,k) = Q.col(k).mean();
        }

        //compute sigma2_s
        MatrixXd sigma2_s(L, 1);
        for (int l = 0; l < L; l++) {
                sigma2_s(l, 0) = (q * (All.row(l).array() * (1 - All.row(l).array())).matrix().transpose()).sum();
        }

        //compute sigma2_T
        MatrixXd sigma2_T(L, 1);
        for (int l = 0; l < L; l++) {
                double aux = (q * All.row(l).transpose()).sum();
                sigma2_T(l, 0) = aux * (1 -aux);
        }

        MatrixXd Fst(L, 1);
        for (int l = 0; l<L; l++) {
          if (sigma2_T(l,0) != 0.0) {
            Fst(l,0) = 1 - sigma2_s(l,0) / sigma2_T(l,0);
          } else {
            //to avoid numerical issues
            Fst(l,0) = 1;
          }
        }

        //Compute Fst
        return (Fst);
}
