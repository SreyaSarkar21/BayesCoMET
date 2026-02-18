#' @title comet
#'
#' @description This function implements the collapsed Gibbs sampler for fitting the compressed mixed-effects tensor (CoMET) model.
#'
#' @param yijs a vector containing the response for all the observations.
#' @param xijlist a list, each component contains a tensor-valued fixed-effect covariate for each observation.
#' @param zijlist a list, each component contains a tensor-valued random-effect covariate for each observation.
#' @param mis a vector of cluster sizes.
#' @param K a user-specified rank for CP structure of fixed-effect coefficient. Should not exceed the dimensions for matrix-valued covariates.
#' @param kdims a vector of length \eqn{D} (number of tensor modes), where each element equals the \eqn{d}-th mode-specific compressed covariance dimension \eqn{k_d}.
#' @param a0 shape hyperparameter for inverse-gamma prior for idiosyncratic error variance.
#' @param b0 scale hyperparameter for inverse-gamma prior for idiosyncratic error variance.
#' @param gammaVar0 a vector of length \eqn{D} (number of tensor modes) with each element being the variance hyperparameter for the Gaussian prior on the columnwise vectorization of the mode-specific compressed covariance parameter.
#' @param R_list a list of length \eqn{D} (number of tensor modes). Each component corresponds to a \eqn{k_d \times q} random projection matrix (for compressing the random slope) with iid entries following \eqn{\mathcal{N}(0, 1/k_d)}, where \eqn{k_d} is the compressed covariance dimension for mode-\eqn{d}.
#' @param S_list a list of length \eqn{D} (number of tensor modes). Each component corresponds to a \eqn{k_d \times q} random projection matrix (for compressing the random-effect covariate) with iid entries following \eqn{\mathcal{N}(0, 1/k_d)}, where \eqn{k_d} is the compressed covariance dimension for mode-\eqn{d}.
#' @param niter number of MCMC iterations.
#' @param nburn number of burn-in samples.
#' @param nthin thinning size specifying every nthin-th sample is retained.
#' @return a list containing the following components:
#' \describe{
#' \item{betaSamp}{a matrix with \eqn{\bigl(\text{niter}-\text{nburn}\bigr)/\text{nthin}} rows and \eqn{p} columns containing posterior samples of fixed effect parameter \eqn{\beta}.}
#' \item{errVarSamp}{a numeric vector of length \eqn{\bigl(\text{niter}-\text{nburn}\bigr)/\text{nthin}} containing posterior samples of error variance \eqn{\tau^2}.}
#' \item{gammaSamp}{a matrix with \eqn{\bigl(\text{niter}-\text{nburn}\bigr)/\text{nthin}} rows and \eqn{k_1 k_2} columns containing posterior samples of columnwise vectorization of the compressed covariance parameter \eqn{\Gamma}.}
#' \item{sampler_time}{time taken by the sampler to run.}
#' }
#' @import Matrix
#' @import rTensor
#' @importFrom stats rnorm rgamma quantile
#' @export


comet <- function(yijs, xijlist, zijlist, mis, K, kdims,
                  a0, b0, gammaVar0, R_list, S_list,
                  niter, nburn, nthin) {

    n <- length(mis) ## number of groups/clusters
    N <- sum(mis) ## total number of observations
    pdims <- dim(xijlist[[1]]) ## fixed-effect tensor dimensions
    qdims <- dim(zijlist[[1]]) ## random-effect tensor dimensions
    nmodes <- length(pdims) ## number of tensor modes

    stopifnot(length(yijs) == N)
    stopifnot(length(xijlist) == N)
    stopifnot(length(zijlist) == N)
    if(nmodes == 2 & K > min(pdims)) {
        stop("`K` must be less than or equal to matrix dimensions.")
    }
    stopifnot(length(kdims) == nmodes)
    stopifnot(length(gammaVar0) == nmodes)
    stopifnot(length(R_list) == nmodes)
    stopifnot(length(S_list) == nmodes)
    stopifnot(a0 > 0); stopifnot(b0 > 0)

    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)

    ### Prior covariance matrices for compressed covariance parameters ###
    Sigma_gamma_list <- lapply(1:nmodes, function(d) {gammaVar0[d] * diag(kdims[d] * kdims[d])})

    ################## Storing final samples #####################
    betaSamp <- matrix(NA, (niter - nburn) / nthin, prod(pdims)) ## vectorized fixed-effect coefficient
    errVarSamp <- rep(NA, (niter - nburn) / nthin) ## idiosyncratic error variance
    lambda2Samplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)}) ## local shrinkage parameters
    delta2Samp <- matrix(NA, (niter - nburn) / nthin, K) ## global shrinkage parameters
    nuSamplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)})
    xiSamp <- matrix(NA, (niter - nburn) / nthin, K)
    gammaSamplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, kdims[d] * kdims[d])}) ## vectorized compressed covariance parameters
    #############################################################

    ####################### Initialization ######################
    Gamma_list <- lapply(1:nmodes, function(d) {matrix(rnorm(kdims[d] * kdims[d], 0, sd = sqrt(gammaVar0[d])), kdims[d], kdims[d])})
    errVar <- 1 / rgamma(1, shape = a0, rate = b0)
    nu_list <- lapply(1:nmodes, function(d) {matrix(1/rgamma(pdims[d] * K, shape = 0.5, rate = 1), pdims[d], K)})
    lambda2_list <- lapply(1:nmodes, function(d) {matrix(1/rgamma(pdims[d] * K, shape = 0.5, rate = 1/as.vector(nu_list[[d]])), pdims[d], K)})
    xi_list <- 1 / rgamma(K, shape = 0.5, rate = 1)
    delta2_list <- 1 / rgamma(K, shape = 0.5, rate = 1/xi_list)

    ### ridge-type initialization for fixed-effect coefficient B ###
    x_rows <- lapply(xijlist, as.vector) ## Vectorize tensors X_ij's
    x_mat <- do.call(rbind, x_rows)
    gram_mat <- crossprod(x_mat)
    # OLS or ridge-like estimate
    beta_init <- chol2inv(chol(gram_mat + diag(1e-6, nrow(gram_mat)))) %*% crossprod(x_mat, yijs)
    #message("Done beta_init")
    B <- array(beta_init, dim = pdims)
    B_factors <- init_CP_factors(beta_vec = beta_init, pdims = pdims, K = K)$U
    ###########################################################

    ### compute compressed random-effect tensor covariates \mathcal{\tilde Z}_{ij} ###
    comp_zijlist <- vector("list", N)
    Skron_not1 <- revkronLOO(S_list, d = 1)
    for(ij in 1:N) {
        zij_tilde_mode1 <- S_list[[1]] %*% tcrossprod(mode_matricize(zijlist[[ij]], d = 1), Skron_not1)
        comp_zijlist[[ij]] <- array(zij_tilde_mode1, dim = kdims)
    }
    vec_comp_zijlist <- lapply(comp_zijlist,
                               function(foo) {as.vector(foo)})
    zi_tilde_list <- lapply(1:n,
                            function(i) {
                                rows_i <- mis_starts[i]:mis_cumsum[i]
                                do.call("rbind", vec_comp_zijlist[rows_i])
                            })

    vecxijlist <- lapply(xijlist, as.vector)

    RRtkron <- revkronAll(lapply(R_list, tcrossprod))

    cts <- 0
    startTime <- proc.time()
    for(its in 1:niter) {
        if(its %% 1000 == 0) cat("iteration: ", its, "\n")
        cycle1Samp <- cometCycle1(yijs = yijs, vecxijlist = vecxijlist, comp_zijlist = comp_zijlist,
                                   zi_tilde_list = zi_tilde_list, mis = mis,
                                   B = B, errVar = errVar,
                                   Gamma_list = Gamma_list, Sigma_gamma_list = Sigma_gamma_list,
                                   kdims = kdims, RRtkron = RRtkron)
        ### Update compressed covariance parameters Gammas ###
        Gamma_list <- cycle1Samp$GammaSamplist

        Gkron <- revkronAll(Gamma_list)
        cycle2Samp <- cometCycle2(yijs = yijs, xijlist = xijlist,
                                   zi_tilde_list = zi_tilde_list,
                                   mis = mis,
                                   B_factors = B_factors, K = K,
                                   a0 = a0, b0 = b0, errVar = errVar,
                                   lambda2_list = lambda2_list, delta2_list = delta2_list,
                                   nu_list = nu_list, xi_list = xi_list,
                                   Gkron = Gkron, RRtkron = RRtkron)

        ###### Update all other parameters ######
        for(d in 1:nmodes) {
            B_factors[[d]] <- cycle2Samp$BfactorsSamp[[d]]
            lambda2_list[[d]] <- cycle2Samp$lambda2Samplist[[d]]
            nu_list[[d]] <- cycle2Samp$nuSamplist[[d]]
        }
        B <- cycle2Samp$BSamp
        xi_list <- cycle2Samp$xiSamplist
        delta2_list <- cycle2Samp$delta2Samplist
        errVar <- cycle2Samp$errVarSamp

        ### Store samples ###
        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            for(d in 1:nmodes) {
                gammaSamplist[[d]][cts, ] <- as.vector(cycle1Samp$GammaSamplist[[d]])
                lambda2Samplist[[d]][cts, ] <- as.vector(cycle2Samp$lambda2Samplist[[d]])
                nuSamplist[[d]][cts, ] <- as.vector(cycle2Samp$nuSamplist[[d]])
            }
            betaSamp[cts, ] <- as.vector(cycle2Samp$BSamp)
            errVarSamp[cts] <- cycle2Samp$errVarSamp
            delta2Samp[cts, ] <- cycle2Samp$delta2Samplist
            xiSamp[cts, ] <- cycle2Samp$xiSamp
        }
    }
    endTime <- proc.time()

    list(betaSamp = betaSamp, errVarSamp = errVarSamp,
         gammaSamplist = gammaSamplist,
         sampler_time = endTime - startTime)
}
