#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Returns the mode-n unfolding as an arma::mat.
// A is a NumericVector with dim attribute (like an R array).
// d denotes the mode index to be matricized.
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

// Returns G_n kronecker ... G_1 as an arma::mat for a list of n matrices.
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


// Returns G_n kronecker G_{d+1} kronecker G_{d-1} ... G_1 as an arma::mat
// Glist is a list of n matrices and d is the mode index to be left out.
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

// matching-columnwise Khatri-Rao product of two matrices with same number of columns
// [[Rcpp::export]]
arma::mat khatri_rao_comet_cpp(const arma::mat& A, const arma::mat& B) {
    if (A.n_cols != B.n_cols) stop("khatri_rao_comet_cpp: ncol(A) != ncol(B)");
    const arma::uword I = A.n_rows, J = B.n_rows, K = A.n_cols;
    arma::mat out(I * J, K, arma::fill::none);
    
    for(arma::uword k = 0; k < K; ++k) {
        // columnwise kronecker: kron(A[, k], B[, k])
        arma::uword idx = 0;
        for(arma::uword j = 0; j < J; ++j) {
            const double bj = B(j, k);
            for(arma::uword i = 0; i < I; ++i) {
                out(idx++, k) = A(i, k) * bj;
            }
        }
    }
    return out;
}

// CP reconstruction
// factorsList: list of factor matrcies, each of dim pd x K
// returns a NumericVector with dim attribute (tensor)
// [[Rcpp::export]]
NumericVector B_cp_cpp(const List& factorsList) {
    const int N = factorsList.size();
    if (N < 2) stop("B_cp_cpp: Need >= 2 factor matrices.");
    
    arma::mat M0 = as<arma::mat>(factorsList[0]);
    const arma::uword K = M0.n_cols;
    
    IntegerVector dims(N);
    dims[0] = (int)M0.n_rows;
    for (int d = 1; d < N; ++d) {
        arma::mat Md = as<arma::mat>(factorsList[d]);
        if (Md.n_cols != K) stop("B_cp_cpp: all factors must have same number of columns (K).");
        dims[d] = (int)Md.n_rows;
    }
    
    // Total size
    long long total_ll = 1;
    for (int d = 0; d < N; ++d) total_ll *= (long long)dims[d];
    if (total_ll > std::numeric_limits<R_xlen_t>::max()) stop("B_cp_cpp: tensor too large.");
    
    NumericVector out((R_xlen_t)total_ll);
    out.attr("dim") = dims;
    
    // R column-major: index = i1 + d1*(i2 + d2*(i3 + ...))
    // Accumulate rank-1 components.
    std::vector<long long> stride(N, 1);
    for (int d = 1; d < N; ++d) stride[d] = stride[d - 1] * (long long)dims[d - 1];
    
    std::vector<arma::mat> U(N);
    for (int d = 0; d < N; ++d) U[d] = as<arma::mat>(factorsList[d]);
    
    // Iterate over all entries
    // Compute sum_{r=1..K} prod_{d=1..N} U_d[i_d, r]
    for (R_xlen_t lin = 0; lin < out.size(); ++lin) {
        // decode lin -> multi-index
        long long tmp = (long long)lin;
        std::vector<int> idx(N, 0);
        for (int d = 0; d < N; ++d) {
            idx[d] = (int)(tmp % dims[d]);
            tmp /= dims[d];
        }
        
        double val = 0.0;
        for (arma::uword r = 0; r < K; ++r) {
            double prod = 1.0;
            for (int d = 0; d < N; ++d) prod *= U[d](idx[d], r);
            val += prod;
        }
        out[lin] = val;
    }
    
    return out;
}


