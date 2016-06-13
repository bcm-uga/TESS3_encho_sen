#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]



//' TODO
//' @export
// [[Rcpp::export]]
Rcpp::List SampleGenoFromGenerativeModelTESS3(const Rcpp::NumericMatrix & Q, const Rcpp::NumericMatrix & G, const Rcpp::NumericMatrix & coord, int ploidy, int openMP_core_num = 1) {
  #ifdef _OPENMP
        omp_set_num_threads(openMP_core_num);
  #endif

        int D = ploidy + 1;
        int L = G.nrow() / D;
        int n = Q.nrow();
        int K = Q.ncol();

        // Because Rcpp not understant with the transpose...
        Map<MatrixXd> eQ(as<Map<MatrixXd> >(Q));
        Map<MatrixXd> eG(as<Map<MatrixXd> >(G));
        MatrixXd ePt = eG * eQ.transpose();
        NumericMatrix X(n, L);

        // auxiliary variables
#ifdef _OPENMP
        ArrayXXi ans(D, omp_get_max_threads());
#else
        ArrayXXi ans(D, 1);
#endif
        ArrayXi allele(D);

        //fill allele
        for (int p = 0; p < D; p++) {
                allele(p) = p;
        }
        #pragma omp parallel for
        for( int i = 0; i < n; i++) {
                for( int j = 0; j < L; j++) {
#ifdef _OPENMP
                        R::rmultinom(1, &ePt(j * D, i), D, &ans(0, omp_get_thread_num()));
                        X(i, j) = (allele * ans.col(omp_get_thread_num())).sum();
#else
                        R::rmultinom(1, &ePt(j * D, i), D, &ans(0, 0) );
                        X(i, j) = (allele * ans.col(0)).sum();
#endif
                }
        }

        return Rcpp::List::create(Rcpp::Named("n") = n,
                                  Rcpp::Named("K") = K,
                                  Rcpp::Named("ploidy") = ploidy,
                                  Rcpp::Named("Q") = Q,
                                  Rcpp::Named("G") = G,
                                  Rcpp::Named("coord") = coord,
                                  Rcpp::Named("X") = X);
}

//' TODO
//'
// [[Rcpp::export]]
Eigen::MatrixXi ComputeZHelper(const Eigen::Map<Eigen::MatrixXd> Q,int n, int L) {
        Eigen::MatrixXi Z(n, L);

        for (int i = 0; i < n; i++) {
                for (int j = 0; j < L; j++) {
                        Z(i, j) = Rcpp::rbinom(1, 1, 1 - Q(i,0))[0] + 1; // because sample is not available...
                }
        }
        return Z;
}

class Offset {
private:
int nrows, ncols, nmats;

public:
Offset( int nrows_, int ncols_, int nmats_) : nrows(nrows_),
        ncols(ncols_), nmats(nmats_){
}

int operator()( int i, int j, int k){
        return i + j * nrows + k * ( nrows * ncols );
}

};

//' TODO
//'
// [[Rcpp::export]]
Eigen::MatrixXi ComputeAdmixtedGeno(const Rcpp::NumericVector & geno, const Eigen::Map<Eigen::MatrixXi> Z, int n, int L) {
        Eigen::MatrixXi adGeno(n, L);
        Offset offset( n, L, 2 );
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < L; j++) {
                        adGeno(i, j) = geno[offset(i, j, Z(i,j) - 1)];
                }
        }
        return adGeno;
}
