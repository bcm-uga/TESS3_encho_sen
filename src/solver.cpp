// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <time.h>

using namespace Eigen;
using namespace Rcpp;

//******************************************************************************
//****************************** Helper functions ******************************

// project Q into the constraint space
void ProjectQ(Eigen::Map<Eigen::MatrixXd> Q) {
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
void ProjectG(Eigen::Map<Eigen::MatrixXd> G, int D) {
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

// [[Rcpp::export]]
void InitOpenMP(int n) {
#ifdef _OPENMP
        if (omp_get_num_procs() < n) {
          Rcpp::Rcout << "openMP.core.num is greater than number of cores available. openMP.core.num is set to " << omp_get_num_procs() << std::endl ;
          n = omp_get_num_procs();
        }
        omp_set_num_threads(n);
        Eigen::initParallel();
#endif
        return;
}

//' solve min || X - Q G^T|| + lambda * tr(Q^T Lapl Q)
// [[Rcpp::export]]
Rcpp::List ComputeMCPASolution(const Eigen::Map<Eigen::MatrixXd> X, int K, const Eigen::Map<Eigen::MatrixXd> Lapl, double lambdaPrim, int D, int maxIteration, double tolerance, Eigen::Map<Eigen::MatrixXd> Q, Eigen::Map<Eigen::MatrixXd> G) {
      // time
      clock_t start, end;
      NumericVector times(maxIteration, NA_REAL);


        // Some const
        const int L = X.cols() / D;
        const int n = X.rows();

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
#ifdef _OPENMP
        // multithreaded OpenMP version of code
        Rcpp::Rcout << "Main loop with " << omp_get_max_threads() << " threads: " << std::endl;
#else
        Rcpp::Rcout << "Main loop :" << std::endl;
#endif
        // variables
        int it = 0;
        bool converg = FALSE;
        LDLT<MatrixXd> lltG;
        while (!converg && it < maxIteration)
        #pragma omp parallel
        {
          start = clock();
                // update G
                #pragma omp master
                {
                  lltG = (Q.transpose() * Q).ldlt();
                }

                #pragma omp barrier
                #pragma omp for
                for (int j = 0; j < G.rows(); j++) {
                  G.row(j) = (lltG.solve(Q.transpose() * X.col(j))).transpose();
                }
                #pragma omp master
                {
                  ProjectG(G, D);

                // update Q
                  RQ = R * Q;
                }
                #pragma omp barrier
                #pragma omp for
                for (int i = 0; i < n; i++) {
                        RQ.row(i) = ((G.transpose() * G + lambda * vps(i) * Ik).ldlt().solve(G.transpose() * RX.row(i).transpose())).transpose();
                }
                #pragma omp master
                {
                Q = R.transpose() * RQ;
                ProjectQ(Q);

                // compute normalized residual error
                errAux = (X - Q * G.transpose()).norm() / X.norm();
                // Rcpp::Rcout << "iteration" << it << "& error : " << err << std::endl; // debug
                Rcpp::Rcout << "---iteration: " << it <<"/" << maxIteration << std::endl;
                // Test the convergence
                converg = (std::abs(errAux - err) < tolerance);
                err = errAux;
                end = clock();
                times[it] = ((double) (end - start)) / CLOCKS_PER_SEC;
                it++;
              }

        }
        return Rcpp::List::create(Rcpp::Named("it") = it,
                          Rcpp::Named("times") = times);
}
