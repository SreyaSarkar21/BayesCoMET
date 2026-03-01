#' @title cometCycle1
#'
#' @description This function implements the collapsed Gibbs sampler for fitting the compressed mixed-effects tensor (CoMET) model.
#'
#' @param y a vector containing the response for all the observations.
#' @param vecxlist a list, each component contains the vectorized fixed-effect covariate for each observation.
#' @param comp_zlist a list, each component contains the compressed random-effect covariate tensor \eqn{\tilde{\mathcal{Z}}} for each observation.
#' @param z_tilde_list a list, each component contains vectorized compressed random-effect covariates, stacked for each cluster.
#' @param mis a vector of cluster sizes.
#' @param B \eqn{\mathcal{B}} sampled in previous iteration.
#' @param errVar \eqn{\tau^2} sampled in previous iteration.
#' @param Gamma_list a list of length \eqn{D}, each component contains \eqn{\Gamma_d} sampled in previous iteration.
#' @param Sigma_gamma_list a list of length \eqn{D}, each component corresponds to the prior covariance matrix for the corresponding \eqn{\gamma_d = \operatorname{vec}(\Gamma_d)}.
#' @param kdims vector of compressed covariance dimensions.
#' @param RRtkron \eqn{R_D R_D^{\top} \otimes \dots \otimes R_1 R_1^{\top}} computed at the start of the model fitting.
#' @return a list containing the following components:
#' \describe{
#' \item{GammaSamplist}{a list of \eqn{D} matrices each of dimension \eqn{k_d \times k_d}, denoting the updated compressed covariance parameters at a given iteration.}
#' \item{Di_tilde}{an array of dimension \eqn{k_1 \times \dots \times k_D}, denoting the imputed compressed random slope tensor at a given iteration. }
#' }
#' @importFrom stats rnorm

cometCycle1 <- function(y, vecxlist, comp_zlist, z_tilde_list, mis,
                         B, errVar,
                         Gamma_list, Sigma_gamma_list, kdims, RRtkron) {

    n <- length(mis); N <- sum(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    nmodes <- length(kdims)

    #### calculate residuals given B
    check_y <- vapply(seq_len(N), function(ij) {y[ij] - sum(vecxlist[[ij]] * as.vector(B))}, 0.0)
    check_ylist <- lapply(seq_len(n), function(i) {rows_i <- mis_starts[i]:mis_cumsum[i]; check_y[rows_i]})

    #### Impute compressed random slopes \tilde D_i ####
    di <- vector("list", n)
    Di_tilde <- vector("list", n)
    Gkron <- revkronAll(Gamma_list)
    for (i in 1:n) {
        Vyd <- errVar * z_tilde_list[[i]] %*% Gkron %*% RRtkron
        VyyInv <- chol2inv(chol(Vyd %*% crossprod(Gkron, t(z_tilde_list[[i]])) + errVar * diag(1, nrow(z_tilde_list[[i]]))))
        VydTvyyInv <- crossprod(Vyd, VyyInv)
        postRanVar <- errVar * RRtkron - drop(VydTvyyInv %*% Vyd)
        postRanMean <- drop(VydTvyyInv %*% check_ylist[[i]])
        di[[i]] <- postRanMean + drop(crossprod(chol(postRanVar + diag(1e-6, nrow(postRanVar))), rnorm(length(postRanMean))))
        Di_tilde[[i]] <- array(di[[i]], dim = kdims)
    }

    #### Sample compressed covariance parameters ####
    for(d in 1:nmodes) {
        check_z_gamma_groups <- vector("list", n)
        for(i in 1:n) {
            rows_i <- mis_starts[i]:mis_cumsum[i]
            check_z_gamma_groups[[i]]  <- do.call("rbind", lapply(comp_zlist[rows_i],
                                                                  function(foo) {as.vector(mode_matricize(foo, d) %*% tcrossprod(revkronLOO(Gamma_list, d), mode_matricize(Di_tilde[[i]], d)))}))
        }
        check_z_gamma <- do.call("rbind", check_z_gamma_groups)

        postVar_gamma <- chol2inv(chol(crossprod(check_z_gamma) / errVar +
                                           chol2inv(chol(Sigma_gamma_list[[d]]))))
        postMean_gamma <- drop(postVar_gamma %*% crossprod(check_z_gamma, check_y) / errVar)

        gammaSamp <- postMean_gamma + drop(crossprod(chol(postVar_gamma + diag(1e-6, nrow(postVar_gamma))), rnorm(length(postMean_gamma))))
        Gamma_list[[d]] <- matrix(gammaSamp, kdims[d], kdims[d])
    }

    list(GammaSamplist = Gamma_list, vecDi_tilde = di)
}
