#' @title CoMET
#'
#' @description This function implements the collapsed Gibbs sampler for fitting the compressed mixed-effects tensor (CoMET) model.
#'
#' @param yijs a vector containing the response for all the observations.
#' @param xlist a list where each component contains a tensor-valued fixed-effect covariate for each observation.
#' @param zlist a list where each component contains a tensor-valued random-effect covariate for each observation.
#' @param a0 shape hyperparameter for inverse-gamma prior for idiosyncratic error variance.
#' @param b0 scale hyperparameter for inverse-gamma prior for idiosyncratic error variance.
#' @param gammaVar0 a vector of length equal to the number of tensor modes, each element in this vector corresponds is the variance hyperparameter for the Gaussian prior on the columnwise vectorization of the mode-specific compressed covariance parameter.
#' @param R_list a list of length equal to the number of tensor modes, each component corresponds to a \eqn{k_d \times q} random projection matrix with entries iid from a Gaussian distribution with mean 0 and variance \eqn{1/k_d}, where \eqn{k_d} is the compressed covariance dimension for mode-\eqn{d}.
#' @param S_list a list of length equal to the number of tensor modes, each component corresponds to a \eqn{k_d \times q} random projection matrix with entries iid from a Gaussian distribution with mean 0 and variance \eqn{1/k_d}, where \eqn{k_d} is the compressed covariance dimension for mode-\eqn{d}.
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
 

CoMET <- function(yijs, Xijlist, Zijlist, mis, K,
                          a0, b0, gammaVar0, R_list, S_list,
                          niter, nburn, nthin) {

    n <- length(mis); N <- sum(mis)
    pdims <- dim(Xijlist[[1]]); qdims <- dim(Zijlist[[1]])
    nmodes <- length(pdims)
    kdims <- unlist(lapply(R_list, nrow))
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    
    Sigma_gamma_list <- lapply(1:nmodes, function(d) {gammaVar0[d] * diag(kdims[d] * kdims[d])})
    
    ################## Storing final samples #####################
    gammaSamplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, kdims[d] * kdims[d])})
    # betatildeSamplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)})
    # GkronSamp <- matrix(NA, (niter - nburn) / nthin, prod(kdims) ^ 2)
    betaSamp <- matrix(NA, (niter - nburn) / nthin, prod(pdims))
    errVarSamp <- rep(NA, (niter - nburn) / nthin)
    lambda2Samplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)}) 
    nuSamplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)})
    delta2Samp <- matrix(NA, (niter - nburn) / nthin, K) ## global shrinkage
    xiSamp <- matrix(NA, (niter - nburn) / nthin, K)
    #############################################################
    
    ####################### Initialization ######################
    Gamma_list <- lapply(1:nmodes, function(d) {matrix(rnorm(kdims[d] * kdims[d], 0, sd = sqrt(gammaVar0[d])), kdims[d], kdims[d])})
    #Bfactors <- lapply(1:nmodes, function(d) {matrix(0, pdims[d], K)})
    #B <- B_cp(Bfactors) + array(rep(1e-4, prod(pdims)), dim = pdims)
    #Bfactors <- lapply(1:nmodes, function(d) {matrix(runif(pdims[d] * K, -5, 5), pdims[d], K)})
    #B <- B_cp(Bfactors)
    errVar <- 1 / rgamma(1, shape = a0, rate = b0)
    nu_list <- lapply(1:nmodes, function(d) {matrix(1/rgamma(pdims[d] * K, shape = 0.5, rate = 1), pdims[d], K)})
    lambda2_list <- lapply(1:nmodes, function(d) {matrix(1/rgamma(pdims[d] * K, shape = 0.5, rate = 1/as.vector(nu_list[[d]])), pdims[d], K)})
    xi_list <- 1 / rgamma(K, shape = 0.5, rate = 1)
    delta2_list <- 1 / rgamma(K, shape = 0.5, rate = 1/xi_list)
    
    ######### OLS/Ridge initialization ##########
    x_rows <- lapply(Xijlist, as.vector) ## Vectorize tensors X_ij's
    X_mat <- do.call(rbind, x_rows)
    gram_mat <- crossprod(X_mat)
    # OLS or Ridge-like estimate
    beta_init <- chol2inv(chol(gram_mat + diag(1e-6, nrow(gram_mat)))) %*% crossprod(X_mat, yijs)
    #message("Done beta_init")
    B <- array(beta_init, dim = pdims)
    # if(length(pdims) == 2) {
    #     sv <- svd(B); K <- min(K, min(pdims))
    #     U <- sv$u[, 1:K, drop=FALSE]; V <- sv$v[, 1:K, drop=FALSE]; Dk <- sv$d[1:K]
    #     U <- sweep(U, 2, sqrt(colSums(U^2)), "/")
    #     V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
    #     lambda <- Dk  # after re-normalization, you can absorb norms into lambda as desired
    #     list(U = list(U, V), lambda = lambda)
    # }
    #cp_B <- rTensor::cp(rTensor::as.tensor(B), num_components = K)
    Bfactors <- init_CP_factors(beta_vec = beta_init, pdims = pdims, K = K)$U
    #message("Done B factors init")
    #Bfactors[[1]] <- Bfactors[[1]] %*% diag(cp_B$lambdas, K, K)
    # # HOSVD initialization
    # Bfactors <- list()
    # for (d in 1:length(pdims)) {
    #     unfold <- mode_matrix(B, d)
    #     svd_res <- svd(unfold)
    #     Bfactors[[d]] <- svd_res$u[,1:K] %*% diag(svd_res$d[1:K] ^(1/length(pdims)))
    # }
    ###########################################################
    
    #### compress Zij arrays ####
    comp_Zijlist <- vector("list", N)
    Skron_not1 <- LOORevKron(S_list, d = 1)
    for(ij in 1:N) {
        Zij_tilde_mode1 <- S_list[[1]] %*% tcrossprod(mode_matrix(Zijlist[[ij]], n = 1), Skron_not1)
        comp_Zijlist[[ij]] <- array(Zij_tilde_mode1, dim = kdims)
    }
    #### vectorize the compressed Zij arrays ####
    vec_comp_Zijlist <- lapply(comp_Zijlist,
                               function(foo) {as.vector(foo)})
    ###############################
    Zi_tilde_list <- lapply(1:n,
                            function(i) {
                                rows_i <- mis_starts[i]:mis_cumsum[i]
                                do.call("rbind", vec_comp_Zijlist[rows_i])
                            })
    # #### vectorize the original Zij arrays ####
    # vec_Zijlist <- lapply(Zijlist,
    #                       function(foo) {as.vector(foo)})
    # Stkron <- t(kronAll(S_list))
    # ###############################
    # Zi_tilde_list <- vector("list", n)
    # for(i in 1:n) {
    #     rows_i <- mis_starts[i]:mis_cumsum[i]
    #     Zi_tilde_list[[i]] <- do.call("rbind", vec_Zijlist[rows_i]) %*% Stkron
    # }
    
    cts <- 0
    startTime <- proc.time()
    for(its in 1:niter) {
        if(its %% 1000 == 0) cat("iteration: ", its, "\n")
        cycle1Samp <- cmeCycle1_SG(yijs = yijs, Xijlist = Xijlist, comp_Zijlist = comp_Zijlist,
                                   Zi_tilde_list = Zi_tilde_list, mis = mis,
                                   B = B, errVar = errVar,
                                   Gamma_list = Gamma_list, Sigma_gamma_list = Sigma_gamma_list,
                                   R_list = R_list)
        ###### Update Gamma's ######
        Gamma_list <- cycle1Samp$GammaSamplist
        #print("Cycle1 Done")
        
        cycle2Samp <- cmeCycle2_SG(yijs = yijs, Xijlist = Xijlist,
                                   Zi_tilde_list = Zi_tilde_list,
                                   mis = mis,
                                   B_list = Bfactors, K = K,
                                   a0 = a0, b0 = b0, errVar = errVar,
                                   lambda2_list = lambda2_list, delta2_list = delta2_list,
                                   nu_list = nu_list, xi_list = xi_list,
                                   Gamma_list = Gamma_list, R_list = R_list)
        
        ###### Update all other parameters ######
        for(d in 1:nmodes) {
            Bfactors[[d]] <- cycle2Samp$BfactorsSamp[[d]]
            lambda2_list[[d]] <- cycle2Samp$lambda2Samplist[[d]]
            nu_list[[d]] <- cycle2Samp$nuSamplist[[d]]
        }
        B <- cycle2Samp$BSamp
        xi_list <- cycle2Samp$xiSamplist
        delta2_list <- cycle2Samp$delta2Samplist
        errVar <- cycle2Samp$errVarSamp
        #print("Cycle2 Done")
        
        ### Store samples ###
        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            for(d in 1:nmodes) {
                gammaSamplist[[d]][cts, ] <- as.vector(cycle1Samp$GammaSamplist[[d]])
                # betatildeSamplist[[d]][cts, ] <- as.vector(cycle2Samp$BfactorsSamp[[d]])
                lambda2Samplist[[d]][cts, ] <- as.vector(cycle2Samp$lambda2Samplist[[d]])
                nuSamplist[[d]][cts, ] <- as.vector(cycle2Samp$nuSamplist[[d]])
            }
            #GkronSamp[cts, ] <- as.vector(kronAll(cycle1Samp$GammaSamplist))
            betaSamp[cts, ] <- as.vector(cycle2Samp$BSamp)
            errVarSamp[cts] <- cycle2Samp$errVarSamp
            delta2Samp[cts, ] <- cycle2Samp$delta2Samplist
            xiSamp[cts, ] <- cycle2Samp$xiSamp
        }
    }
    endTime <- proc.time()
    
    list(gammaSamplist = gammaSamplist,
         betaSamp = betaSamp, errVarSamp = errVarSamp,
         lambda2Samplist = lambda2Samplist, nuSamplist = nuSamplist,
         delta2Samp = delta2Samp, xiSamp = xiSamp,
         sampler_time = endTime - startTime)
    
    #GkronSamp = GkronSamp, betaFactorsSamp = betatildeSamplist,
}
