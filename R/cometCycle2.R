#' @title cometCycle2
#'
#' @description This function implements the collapsed Gibbs sampler for fitting the compressed mixed-effects tensor (CoMET) model.
#'
#' @param yijs a vector containing the response for all the observations.
#' @param xijlist a list, each component contains a tensor-valued fixed-effect covariate for each observation.
#' @param zi_tilde_list a list, each component contains vectorized compressed random-effect covariates, stacked for each cluster.
#' @param mis a vector of cluster sizes.
#' @param B_factors fixed-effect factor matrices sampled in previous iteration.
#' @param K a user-specified CP rank of fixed-effect coefficient \eqn{\mathcal{B}}. Should not exceed the dimensions for matrix-valued covariates.
#' @param a0 shape hyperparameter for inverse-gamma prior for idiosyncratic error variance.
#' @param b0 scale hyperparameter for inverse-gamma prior for idiosyncratic error variance.
#' @param errVar idiosyncratic error variance sampled in previous iteration.
#' @param lambda2_list a list of length \eqn{D}, each containing \eqn{K p_d} squared local shrinkage \eqn{\lambda^2} values from previous iteration.
#' @param delta2_list a vector of length \eqn{K}, value of rank-\eqn{g}-specific global shrinkage \eqn{\delta^2} values from previous iteration.
#' @param nu_list a list of length \eqn{D}, each containing \eqn{K p_d} auxiliary variable \eqn{\nu} values from previous iteration.
#' @param xi_list a vector of length \eqn{K}, value of rank-\eqn{g}-specific auxiliary variable \eqn{\xi_g} values from previous iteration.
#' @param Gkron \eqn{\Gamma_D \Gamma_D^{\top} \otimes \dots \otimes \Gamma_1 \Gamma_1^{\top}} computed using the \eqn{\Gamma_d} updates from the previous iteration.
#' @param RRtkron \eqn{R_D R_D^{\top} \otimes \dots \otimes R_1 R_1^{\top}} computed at the start of the model fitting.
#' @return a list containing the following components:
#' \describe{
#' \item{BfactorsSamp}{a list of length \eqn{D}, each component comprising the posterior draw of \eqn{p_d \times K} factor matrix \eqn{B_d}.}
#' \item{BSamp}{an array denoting the updated fixed-effect parameter \eqn{\mathcal{B}}.}
#' \item{errVarSamp}{the updated posterior draw of error variance \eqn{\tau^2}.}
#' \item{lambda2Samplist}{a list of length \eqn{D}, each component containing the updated \eqn{K p_d} squared local shrinkage \eqn{\lambda^2} draws.}
#' \item{delta2Samplist}{a vector of length \eqn{K} storing the updated rank-\eqn{g}-specific global shrinkage \eqn{\delta^2} draws.}
#' \item{nuSamplist}{a list of length \eqn{D}, each component containing the updated \eqn{K p_d} auxiliary variable \eqn{nu} draws.}
#' \item{xiSamplist}{a vector of length \eqn{K} storing the updated rank-\eqn{g}-specific auxiliary variable \eqn{\xi_g} draws.}
#' }
#' @importFrom stats rnorm rgamma

cometCycle2 <- function(yijs, xijlist, zi_tilde_list, mis,
                         B_factors, K, a0, b0,
                         errVar, lambda2_list, delta2_list, nu_list, xi_list,
                         Gkron, RRtkron) {

    n <- length(mis); N <- sum(mis)
    pdims <- dim(xijlist[[1]])
    nmodes <- length(pdims)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)

    ######################## Marginalization ########################
    margVari <- vector("list", n)
    for(i in 1:n) {
        margVari[[i]] <- zi_tilde_list[[i]] %*% Gkron %*% RRtkron %*% crossprod(Gkron, t(zi_tilde_list[[i]])) + diag(1, mis[i])
    }
    margVari_svd <- lapply(margVari, svd)
    margVari_inv_half <- vector("list", n)
    for(i in 1:n) {
        if(mis[i] == 1) {
            margVari_inv_half[[i]] <- (margVari_svd[[i]]$u)*(margVari_svd[[i]]$u)/sqrt(margVari_svd[[i]]$d)
        } else {
            margVari_inv_half[[i]] <- margVari_svd[[i]]$u %*% tcrossprod(diag(1 / sqrt(margVari_svd[[i]]$d)), margVari_svd[[i]]$u)
        }
    }
    ylist <- lapply(1:n,
                    function(i) {
                        rows_i <- mis_starts[i]:mis_cumsum[i]
                        margVari_inv_half[[i]] %*% yijs[rows_i]
                    })
    yvec <- as.matrix(unlist(ylist))
    ########################################################################

    ################# Given \tau^2, sample B in 3 steps #################
    for(d in 1:nmodes) {
        xmat_rowlist <- lapply(xijlist,
                               function(foo) {as.vector(mode_matricize(foo, d) %*% khaRaoLOO_B(B_factors, d))})
        xmat_list <- lapply(1:n,
                            function(i) {
                                rows_i <- mis_starts[i]:mis_cumsum[i]
                                temp <- do.call("rbind", xmat_rowlist[rows_i])
                                margVari_inv_half[[i]] %*% temp
                            })

        xmat <- do.call(rbind, xmat_list) ## design matrix for sampling betatilde_d; d = 1, 2, 3
        priorvar_beta <- as.matrix(Matrix::bdiag(lapply(1:K, function(g) {diag(1 / (delta2_list[g] * lambda2_list[[d]][, g]))} )))
        vinv_beta_post <- crossprod(xmat) + priorvar_beta + diag(1e-6, nrow(priorvar_beta))
        covBeta <- errVar * chol2inv(chol(vinv_beta_post))
        xty <- crossprod(xmat, yvec)
        muBeta <- drop(covBeta %*% xty) / errVar
        betatilde <- muBeta + drop(crossprod(chol(covBeta), rnorm(pdims[d] * K)))
        B_factors[[d]] <- matrix(betatilde, pdims[d], K)
        #message("Done B", d)

        for(g in 1:K) {
            lambda2_list[[d]][, g] <- 1 / rgamma(pdims[d], shape = 1,
                                                 rate = 1 / nu_list[[d]][, g] + (betatilde[((g-1)*pdims[d] + 1):(g*pdims[d])] ^ 2) / (2 * delta2_list[g] * errVar))
            nu_list[[d]][, g] <- 1 / rgamma(pdims[d], shape = 1, rate = 1 + 1 / (lambda2_list[[d]][, g]))
        }
    }

    B <- cpB(B_factors)

    ### Given updated B, \lambda's, sample \delta^2_g and \xi_g ###
    g_sums <- sapply(1:K, function(g) {
        sum(sapply(1:nmodes, function(d) {
            sum(B_factors[[d]][, g]^2 / lambda2_list[[d]][, g])
        }))
    })
    delta2_list <- 1 / rgamma(K, shape = 0.5 * (1 + sum(pdims)),
                              rate = 1 / xi_list + g_sums / (2 * errVar))
    xi_list <- 1 / rgamma(K, shape = 1, rate = 1 + 1 / delta2_list)

    ################# Given updated B, sample \tau^2 #################
    xijB_list <- lapply(xijlist,
                        function(foo) {sum(as.vector(foo) * as.vector(B))})
    mu_list <- lapply(1:n,
                      function(i) {
                          rows_i <- mis_starts[i]:mis_cumsum[i]
                          temp2 <- do.call("rbind", xijB_list[rows_i])
                          margVari_inv_half[[i]] %*% temp2
                      })
    prior_contrib <- 0
    for(d in seq_len(nmodes)) {
        for(g in seq_len(K)) {
            prior_contrib <- prior_contrib + sum(B_factors[[d]][, g]^2 / (lambda2_list[[d]][, g] * delta2_list[g]))
        }
    }
    resids <- yvec - drop(do.call(rbind, mu_list))
    shp_err <- a0 + 0.5 * (N + K * sum(pdims)) ## a0 = prior shape of errVar
    scl_err <- b0 + 0.5 * (sum(resids^2) + prior_contrib) ## b0 = prior scale of errVar
    errVar <- 1 / rgamma(1, shape = shp_err, rate = scl_err)
    #print("errVar Done")
    ####################################################################

    ##################### return samples #####################
    list(BfactorsSamp = B_factors,
         BSamp = B, errVarSamp = errVar,
         lambda2Samplist = lambda2_list, delta2Samplist = delta2_list,
         nuSamplist = nu_list, xiSamplist = xi_list)
}

