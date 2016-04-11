// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>


//******************************************************************************
//****************************** Helper functions ******************************

// project Q into the constraint space
void ProjectQ(arma::mat & Q) {
        double rowSum = 0.0;
        for (int i = 0; i < Q.n_rows; i++) {
                rowSum = 0.0;
                for (int j = 0; j < Q.n_cols; j++) {
                        Q(i,j) = std::max(0.0, Q(i,j));
                        rowSum += Q(i,j);
                }
                if (rowSum != 0.0) {
                        Q.row(i) /= rowSum;
                }
        }
}
// project G into the constraint space
void ProjectG(arma::mat & G, int D) {
        int L = G.n_rows / D;
        double sum = 0.0;
        for (int k = 0; k < G.n_cols; k++) {
                for (int l = 0; l < L; l++) {
                        sum = 0.0;
                        for (int j = 0; j < D; j++) {
                                G(D * l + j, k) = std::max(0.0, G(D * l + j, k));
                                sum += G(D * l + j, k);
                        }
                        if (sum != 0.0) {
                                G.submat(D * l, k, D * (l + 1) - 1, k) /= sum;
                        }
                }
        }
}
//******************************************************************************
//******************************************************************************



// [[Rcpp::depends(RcppArmadillo)]]

//' solve min || X - Q G^T|| + lambda * tr(Q^T Lapl Q)
// [[Rcpp::export]]
List ComputeMCPASolution(const arma::mat & X, int K, const arma::mat & Lapl, double lambdaPrim, int D, int maxIteration, double tolerance) {
        // Some const
        const int L = X.n_cols / D;
        const int n = X.n_rows;

        // Init Q and G
        arma::mat G = arma::zeros<mat>(X.n_cols, K);
        arma::mat Q = arma::randu<mat>(X.n_rows, K);
        ProjectQ(Q);

        // Compute Lapl diagonalization
        Rcpp::Rcout << "Computing spectral decomposition of graph laplacian matrix";
        arma::vec vps(n);
        arma::mat R(n, n);
        arma::eig_sym(vps, R, Lapl)
        R = T.t();
        arma::mat RX = R * X;
        Rcpp::Rcout << ": done" << std::endl;

        // Compute lambda
        double vpMax = vps.max();
        double lambda = 0.0;
        if (vpMax != 0.0) {
          lambda = lambdaPrim * (D * L * n) / (K * n * vpMax);
        }

        // constant
        arma::mat Ik = arma::eye<mat>(K,K);

        // auxiliary variables
        arma::mat RQ = Q;
        double err = -10.0;
        double errAux = 0.0;

        // algo
        int it = 0;
        bool converg = FALSE;
        Rcpp::Rcout << "Main loop: " << std::endl;
        while (!converg && it < maxIteration) {
                // update G
                G = arma::solve((Q.t() * Q),Q.t() * X).t();
                ProjectG(G, D);

                // update Q
                RQ = R * Q;
                for (int i = 0; i < n; i++) {
                        RQ.row(i) = ((G.transpose() * G + lambda * vps(i) * Ik).ldlt().solve(G.transpose() * RX.row(i).transpose())).transpose();
                }
                Q = R.t() * RQ;
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
