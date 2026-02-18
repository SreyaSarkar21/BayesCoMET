#' @title mode_matrix
#' @description mode-n matricization of a tensor
#' @param A a numeric array with `dim` attribute.
#' @param d integer, the mode index to be matricized.
#' @return a matrix corresponding to the mode-\eqn{d} unfolding.
#' @export
mode_matrix <- function(A, d) {
    mode_matrix_cpp(A, d)
}

#' @title revkronAll
#' @description Returns \eqn{G_n \otimes \dots G_1} for a list of n matrices.
#' @param Glist a list of matrices
#' @return a matrix yielding kronecker product in reverse order.
#' @export
revkronAll <- function(Glist) {
    revkronAll_cpp(Glist)
}

#' @title LOOrevkron
#' @description Returns \eqn{G_n \otimes G_{d+1} \otimes G_{d-1} \dots G_1}
#' @param Glist a list of matrices
#' @param d an integer denoting the mode to be left out.
#' @return a matrix yielding kronecker product in reverse order.
#' @export
LOOrevkron <- function(Glist, d) {
    LOOrevkron_cpp(Glist, d)
}

#' @title khatri_rao_comet
#' @description matching-columnwise Khatri-Rao product
#' @param A a matrix of dimensions \eqn{p_1 \times k}.
#' @param B a matrix of dimensions \eqn{p_2 \times k}.
#' @return a matrix of dimensions \eqn{p_1 p_2 \times k}.
#' @export
khatri_rao_comet <- function(A, B) {
    khatri_rao_comet_cpp(A, B)
}

B_cp <- function(factorsList) {
    B_cp_cpp(factorsList)
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
        result <- khatri_rao_comet(result, B_list[[i]])
    }
    result
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

