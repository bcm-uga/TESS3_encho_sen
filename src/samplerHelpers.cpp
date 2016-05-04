#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// Trop compliqué...
// Rcpp::List SampleGenoFromGenerativeModelTESS3(const Eigen::Map<Eigen::MatrixXd> Q, const Eigen::Map<Eigen::MatrixXd> G, const Eigen::Map<Eigen::MatrixXd> coord, int ploidy) {
//         int D = ploidy + 1;
//         int L = G.rows() / D;
//         int n = Q.rows();
//         int K = Q.cols();
//         Eigen::MatrixXd P = Q * G.transpose();
//         Eigen::MatrixXd X = P;
//         Rcpp::IntegerVector ans(D);
//         Rcpp::IntegerVector allele(D);
//         // fill allele
//         for (int p = 0; p < D; p++) {
//           allele(p) = p;
//         }
//         for( int i = 0; i < n; i++) {
//                 for( int j = 0; j < L; j++) {
//                     ans = 0;
//                     R::rmultinom(1, P.block(i, j * D, 1, D).data(), D, ans.begin());
//                         X(i, j) = Rcpp::sum(allele * ans);
//                 }
//         }
//
//         return Rcpp::List::create(Rcpp::Named("n") =n,
//                                   Rcpp::Named("K") = K,
//                                   Rcpp::Named("ploidy") = ploidy,
//                                   Rcpp::Named("Q") = Q,
//                                   Rcpp::Named("G") = G,
//                                   Rcpp::Named("coord") = coord,
//                                   Rcpp::Named("X") = X);
// }

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
