#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// --------- mode_matrix (general N-way) ---------
// Returns the mode-n unfolding as an arma::mat.
// A is a NumericVector with dim attribute (like an R array).
// n is 1-based mode index.
// [[Rcpp::export]]
arma::mat mode_matrix_cpp(const NumericVector& A, const int d) {
    IntegerVector dims = A.attr("dim");
    const int n = dims.size();
    if(n < 2) stop("mode_matrix_cpp: input must be at least a 2-dimensional array.");
    if(d < 1 || d > n) stop("mode_matrix_cpp: invalid d.");
    
    // R-arrays are column major
    const int mode = d - 1;
    int nrow = dims[mode];
    long long ncol_ll = 1;
    for(int k = 0; k < n; k++) {
        if(k != mode) {
            ncol_ll *= (long long)dims[k];
        }
    }
    if (ncol_ll > std::numeric_limits<int>::max()) stop("mode_matrix_cpp: too large.");
    int ncol = (int)ncol_ll;
    arma::mat out(nrow, ncol, arma::fill::none);
    
    // Precompute strides for R column-major indexing
    std::vector<long long> stride(n, 1);
    for (int k = 1; k < n; ++k) stride[k] = stride[k - 1] * (long long)dims[k - 1];
    
    // For each linear index in the "other modes" (column of unfolding),
    // map to multi-index excluding 'mode', then fill rows by varying mode index.
    // Build mapping from col -> indices of other dims.
    std::vector<int> other_dims;
    other_dims.reserve(n - 1);
    for (int k = 0; k < n; ++k) if (k != mode) other_dims.push_back(k);
    
    for (int col = 0; col < ncol; ++col) {
        // decode col into indices for other dims
        long long tmp = col;
        std::vector<int> idx(n, 0);
        for (int j = 0; j < (int)other_dims.size(); ++j) {
            int dimj = dims[other_dims[j]];
            idx[other_dims[j]] = (int)(tmp % dimj);
            tmp /= dimj;
        }
        // finally vary mode index for rows
        for (int r = 0; r < nrow; ++r) {
            idx[mode] = r;
            long long lin = 0;
            for (int k = 0; k < n; ++k) lin += (long long)idx[k] * stride[k];
            out(r, col) = A[(R_xlen_t)lin];
        }
    }
    return out;
}


// [[Rcpp::export]]
arma::mat revkronAll_cpp(const List& Glist) {
    const int n = Glist.size();
    if(n < 2) stop("revkronAll_cpp: at least 2 matrices are needed.");
    
    arma::mat result = as<arma::mat>(Glist[n - 1]);
    for(int i = n - 2; i >=0; --i) {
        arma::mat Gi = as<arma::mat>(Glist[i]);
        result = arma::kron(result, Gi);
    }
    return result;
}


// [[Rcpp::export]]
arma::mat LOOrevkron_cpp(const List& Glist, const int d) {
    const int n = Glist.size();
    if(d < 1 || d > n) stop("LOOrevkron_cpp: invalid d.");
    
    arma::mat out(1, 1, arma::fill::ones); // scalar 1
    for(int i = n; i >= 1; --i) {
        if(i == d) continue;
        arma::mat Gi = as<arma::mat>(Glist[i - 1]);
        out = arma::kron(out, Gi);
    }
    return out;
}

