// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

//******************************************************************************
//****************************** Helper functions ******************************

// project Q into the constraint space
void ProjectQ(Eigen::MatrixXd & Q) {
        double rowSum = 0.0;
        for (int i = 0; i < Q.rows(); i++) {
                rowSum = 0.0;
                for (int j = 0; j < Q.cols(); j++) {
                        Q(i,j) = std::max(0.0, Q(i,j));
                        rowSum += Q(i,j);
                }
                if (rowSum != 0.0) {
                        Q.row(i) /= rowSum;
                }
        }
}
// project G into the constraint space
void ProjectG(Eigen::MatrixXd & G, int D) {
        int L = G.rows() / D;
        double sum = 0.0;
        for (int k = 0; k < G.cols(); k++) {
                for (int l = 0; l < L; l++) {
                        sum = 0.0;
                        for (int j = 0; j < D; j++) {
                                G(D * l + j, k) = std::max(0.0, G(D * l + j, k));
                                sum += G(D * l + j, k);
                        }
                        if (sum != 0.0) {
                                G.block(D * l,k,D,1) /= sum;
                        }
                }
        }
}
//******************************************************************************
//******************************************************************************



// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                 // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision

//' solve min || X - Q G^T|| + lambda * tr(Q^T Lapl Q)
// [[Rcpp::export]]
List ComputeMCPASolution(const Eigen::Map<Eigen::MatrixXd> X, int K, const Eigen::Map<Eigen::MatrixXd> Lapl, double lambdaPrim, int D, int maxIteration, double tolerance) {
        // Some const
        const int L = X.cols() / D;
        const int n = X.rows();

        // Init Q and G
        Eigen::MatrixXd G = MatrixXd::Zero(X.cols(), K);
        Eigen::MatrixXd Q = MatrixXd::Random(X.rows(), K);
        Q = Q.cwiseAbs();
        ProjectQ(Q);

        // Compute Lapl diagonalization
        Rcpp::Rcout << "Computing spectral decomposition of graph laplacian matrix";
        SelfAdjointEigenSolver<MatrixXd> es(Lapl);
        VectorXd vps = es.eigenvalues();
        MatrixXd R = es.eigenvectors().transpose();
        MatrixXd RX = R * X;
        Rcpp::Rcout << ": done" << std::endl;

        // Compute lambda
        double vpMax = vps.maxCoeff();
        double lambda = 0.0;
        if (vpMax != 0.0) {
          lambda = lambdaPrim * (D * L * n) / (K * n * vpMax);
        }

        // constant
        MatrixXd Ik = MatrixXd::Identity(K,K);

        // auxiliary variables
        MatrixXd RQ = Q;
        double err = -10.0;
        double errAux = 0.0;

        // algo
        int it = 0;
        bool converg = FALSE;
        Rcpp::Rcout << "Main loop: " << std::endl;
        while (!converg && it < maxIteration) {
                // update G
                G = ((Q.transpose() * Q).ldlt().solve(Q.transpose() * X)).transpose();
                ProjectG(G, D);

                // update Q
                RQ = R * Q;
                for (int i = 0; i < n; i++) {
                        RQ.row(i) = ((G.transpose() * G + lambda * vps(i) * Ik).ldlt().solve(G.transpose() * RX.row(i).transpose())).transpose();
                }
                Q = R.transpose() * RQ;
                ProjectQ(Q);

                // compute normalized residual error
                errAux = (X - Q * G.transpose()).norm() / X.norm();
                // Rcpp::Rcout << "iteration" << it << "& error : " << err << std::endl; // debug
                Rcpp::Rcout << "---iteration: " << it <<"/" << maxIteration << std::endl;
                // Test the convergence
                converg = (std::abs(errAux - err) < tolerance);
                err = errAux;
                it++;
        }


        return List::create(Named("Q") = Q,
                            Named("G") = G);

}
