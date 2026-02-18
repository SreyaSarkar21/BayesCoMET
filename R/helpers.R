#' @title mode_matricize
#' @description mode-n matricization of a tensor
#' @param A a numeric array with `dim` attribute.
#' @param d integer, the mode index to be matricized.
#' @return a matrix corresponding to the mode-\eqn{d} unfolding.
#' @export
mode_matricize <- function(A, d) {
    mode_matricize_cpp(A, d)
}


#' @title revkronAll
#' @description Given \eqn{n} matrices \eqn{G_1, \dots, G_n}, computes \eqn{G_n \otimes \dots \otimes G_1}.
#' @param Glist a list of matrices.
#' @return a matrix yielding Kronecker product of in reverse order.
#' @export
revkronAll <- function(Glist) {
    revkronAll_cpp(Glist)
}


#' @title revkronLOO
#' @description Given \eqn{n} matrices \eqn{G_1, \dots, G_n}, computes \eqn{G_n \otimes \dots \otimes G_{d+1} \otimes G_{d-1} \otimes \dots \otimes G_1}.
#' @param Glist a list of matrices.
#' @param d integer, index corresponding to the matrix to be excluded from the Kronecker product.
#' @return a matrix yielding Kronecker product in reverse order, leaving the \eqn{d}-th matrix.
#' @export
revkronLOO <- function(Glist, d) {
    revkronLOO_cpp(Glist, d)
}


#' @title khatri_rao_comet
#' @description matching-columnwise Khatri-Rao product of two matrices with same number of columns.
#' @param A a \eqn{p_1 \times k} matrix.
#' @param B a \eqn{p_2 \times k} matrix.
#' @return a \eqn{p_1 p_2 \times k} matrix with the \eqn{d}-th column as Kronecker product of the \eqn{d}-th columns of \eqn{A} and \eqn{B}.
#' @export
khatri_rao_comet <- function(A, B) {
    khatri_rao_comet_cpp(A, B)
}


#' @title cpB
#' @description Given a list of factor matrices, each of dimension \eqn{p_d \times K}, reconstructs a \eqn{p_1 \times \dots \times p_D} tensor \eqn{\mathcal{B}} using rank-\eqn{K} CP/PARAFAC structure.
#' @param factorsList list of factor matrices.
#' @return a tensor
#' @export
cpB <- function(factorsList) {
    cpB_cpp(factorsList)
}


#### Computes Khatri-Rao product of the factor matrices of B except mode-d ####
khRaoBLOO <- function(B_list, d) {
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

