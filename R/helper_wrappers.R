mode_matrix <- function(A, d) {
    mode_matrix_cpp(A, d)
}

#### matching-columnwise Khatri-Rao product ####
khatri_rao_cme <- function(A, B) {
    stopifnot(ncol(A) == ncol(B))
    K = matrix(NA, nrow(A) * nrow(B), ncol(A))
    for(ii in 1:ncol(K)) {
        K[, ii] <- kronecker(A[, ii], B[, ii])
    }
    K
}

B_cp <- function(factorsList) {
    N <- length(factorsList)  # Number of modes
    dims <- sapply(factorsList, nrow)
    K <- ncol(factorsList[[1]])
    
    # Check that all have the same number of columns
    if (any(sapply(factorsList, ncol) != K)) {
        stop("All factor matrices must have the same number of columns (i.e., rank K).")
    }
    
    # Initialize output tensor
    B <- array(0, dim = dims)
    
    for (r in 1:K) {
        # Extract the r-th column vectors from each mode factor matrix
        vecs <- lapply(factorsList, function(M) M[, r])
        
        # Build rank-1 component via successive outer products
        rank1 <- vecs[[1]]
        for (i in 2:N) {
            rank1 <- outer(rank1, vecs[[i]])
        }
        
        # Reshape and accumulate
        B <- B + array(rank1, dim = dims)
    }
    B
}

# #### Compute Khatri-Rao product of the factor matrices of B except dim-d ####
B_lookr <- function(B_list, d) {
    if (!is.list(B_list) || length(B_list) < 2) {
        stop("B_list must be a list of at least 2 matrices.")
    }
    
    R_vals <- sapply(B_list, ncol)
    if (length(unique(R_vals)) != 1) {
        stop("All matrices in B_list must have the same number of columns.")
    }
    
    # Get indices excluding d, in reverse order
    rev_indices <- rev(setdiff(seq_along(B_list), d))
    
    # Start with the last matrix in the reversed order
    result <- B_list[[rev_indices[1]]]
    
    # Iterate over the remaining indices in reverse order
    for (i in rev_indices[-1]) {
        result <- khatri_rao_cme(result, B_list[[i]])
    }
    result
}


revkronAll <- function(Glist) {
    revkronAll_cpp(Glist)
}

LOOrevkron <- function(Glist, d) {
    LOOrevkron_cpp(Glist, d)
}


init_CP_factors <- function(beta_vec, pdims, K) {
    B <- array(beta_vec, dim = pdims)
    if (length(pdims) == 2 || any(pdims == 1)) {
        sv <- svd(drop(B)) #; K <- min(K, min(pdims))
        U <- sv$u[, 1:K, drop = FALSE]   # orthonormal
        V <- sv$v[, 1:K, drop = FALSE]   # orthonormal
        lambda <- sv$d[1:K]
        U <- U %*% diag(lambda, nrow = K)
        list(U = list(U, V), lambda = lambda)
    } else {
        cp_B <- rTensor::cp(rTensor::as.tensor(B), num_components = K)
        Bfactors <- cp_B$U
        Bfactors[[1]] <- Bfactors[[1]] %*% diag(cp_B$lambdas, K, K)
        list(U = cp_B$U, lambda = cp_B$lambdas)
    }
}

