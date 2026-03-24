###### Oracle (true covariance structure known) ######

sampleB_SG <- function(yijs, Xijlist, Zi_tilde_list, mis,
                       B_list, K, a0, b0,
                       errVar, lambda2_list, delta2_list, nu_list, xi_list,
                       Sigma) {
    
    pdims <- dim(Xijlist[[1]])
    nmodes <- length(pdims)
    n <- length(mis); N <- sum(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    
    #### Marginalization ####
    margVari <- vector("list", n)
    for(i in 1:n) {
        margVari[[i]] <- Zi_tilde_list[[i]] %*% tcrossprod(Sigma, Zi_tilde_list[[i]]) + diag(1, nrow(Zi_tilde_list[[i]]))
    }
    margVari_svd <- lapply(margVari, svd)
    margVari_inv_half <- vector("list", n)
    for(i in 1:n) {
        if(nrow(Zi_tilde_list[[i]]) == 1) {
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
    
    #### Given \tau^2, sample for B in 3 steps ####
    ################# Given \tau^2, sample B in 3 steps #################
    for(d in 1:nmodes) {
        xmat_rowlist <- lapply(Xijlist,
                               function(foo) {as.vector(BayesCoMET::mode_matricize(foo, d) %*% BayesCoMET:::khaRaoLOO_B(B_list, d))})
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
        B_list[[d]] <- matrix(betatilde, pdims[d], K)
        
        for(g in 1:K) {
            lambda2_list[[d]][, g] <- 1 / rgamma(pdims[d], shape = 1,
                                                 rate = 1 / nu_list[[d]][, g] + (betatilde[((g-1)*pdims[d] + 1):(g*pdims[d])] ^ 2) / (2 * delta2_list[g] * errVar))
            nu_list[[d]][, g] <- 1 / rgamma(pdims[d], shape = 1, rate = 1 + 1 / (lambda2_list[[d]][, g]))
        }
    }
    
    B <- BayesCoMET::cpB(B_list)
    
    ### Given updated B, \lambda's, sample \delta^2_g and \xi_g ###
    g_sums <- sapply(1:K, function(g) {
        sum(sapply(1:nmodes, function(d) {
            sum(B_list[[d]][, g]^2 / lambda2_list[[d]][, g])
        }))
    })
    delta2_list <- 1 / rgamma(K, shape = 0.5 * (1 + sum(pdims)),
                              rate = 1 / xi_list + g_sums / (2 * errVar))
    xi_list <- 1 / rgamma(K, shape = 1, rate = 1 + 1 / delta2_list)
    
    ################# Given updated B, sample \tau^2 #################
    XijB_list <- lapply(Xijlist,
                        function(foo) {sum(foo * B)})
    mu_list <- lapply(1:n,
                      function(i) {
                          rows_i <- mis_starts[i]:mis_cumsum[i]
                          temp2 <- do.call("rbind", XijB_list[rows_i])
                          margVari_inv_half[[i]] %*% temp2
                      })
    prior_contrib <- 0
    for(d in seq_len(nmodes)) {
        for(g in seq_len(K)) {
            prior_contrib <- prior_contrib + sum(B_list[[d]][, g]^2 / (lambda2_list[[d]][, g] * delta2_list[g]))
        }
    }
    resids <- yvec - drop(do.call(rbind, mu_list))
    shp_err <- a0 + 0.5 * (N + K * sum(pdims)) ## a0 = prior shape of errVar
    scl_err <- b0 + 0.5 * (sum(resids^2) + prior_contrib) ## b0 = prior scale of errVar
    errVar <- 1 / rgamma(1, shape = shp_err, rate = scl_err)
    
    #### return samples ####
    list(BfactorsSamp = B_list,
         BSamp = B, errVarSamp = errVar,
         delta2Samplist = delta2_list, lambda2Samplist = lambda2_list, 
         nuSamplist = nu_list, xiSamplist = xi_list)
}

sample_bi_oracle <- function(yijs, Xijlist, Zi_tilde_list, mis, B, Sigma, errVar,
                             jitter = 1e-8) {
    
    n <- length(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    
    q <- ncol(Zi_tilde_list[[1]])
    stopifnot(nrow(Sigma) == q, ncol(Sigma) == q)
    
    # Precompute mean contributions mu_ij = <X_ij, B>
    Bvec <- as.vector(B)
    XijB <- vapply(Xijlist, function(Xij) sum(as.vector(Xij) * Bvec), numeric(1))
    
    bi_mat <- matrix(NA_real_, nrow = n, ncol = q)
    
    for (i in seq_len(n)) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        
        yi <- yijs[rows_i]
        mui <- XijB[rows_i]
        Zi <- Zi_tilde_list[[i]]  # m_i x q
        mi <- length(rows_i)
        
        # Vyy = (Zi Sigma Zi^T + I_mi)
        # adding jitter for numerical stability
        Vyy <- Zi %*% tcrossprod(Sigma, Zi) + diag(1, mi)
        if (jitter > 0) Vyy <- Vyy + diag(jitter, mi)
        
        VyyInv  <- chol2inv(chol(Vyy))
        resid <- yi - mui
        m_post <- as.vector(tcrossprod(Sigma, Zi) %*% (VyyInv %*% resid))
        
        # compute middle = Zi^T Vyy^{-1} Zi
        mid <- crossprod(Zi, VyyInv %*% Zi)
        V_post <- errVar * (Sigma - Sigma %*% mid %*% Sigma)
        
        # Symmetrize & jitter to ensure numerical PSD/PD for chol
        V_post <- (V_post + t(V_post)) / 2
        if (jitter > 0) V_post <- V_post + diag(jitter, q)
        
        cholV <- chol(V_post)
        bi_mat[i, ] <- m_post + as.vector(t(cholV) %*% rnorm(q))
    }
    
    bi_mat
}

samplerOracle_SG <- function(yijs, Xijlist, Zijlist, mis,
                             K, Sigma_list, a0, b0,
                             niter, nburn, nthin) {
    pdims <- dim(Xijlist[[1]])
    qdims <- dim(Zijlist[[1]])
    nmodes <- length(pdims)
    n <- length(mis); N <- sum(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    
    Sigma <- BayesCoMET::revkronAll(Sigma_list)
    
    ################## Storing final samples #####################
    betaSamp <- matrix(NA, (niter - nburn) / nthin, prod(pdims))
    errVarSamp <- rep(NA, (niter - nburn) / nthin)
    lambda2Samplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)}) 
    nuSamplist <- lapply(1:nmodes, function(d) {matrix(NA, (niter - nburn) / nthin, pdims[d] * K)})
    delta2Samp <- matrix(NA, (niter - nburn) / nthin, K)
    xiSamp <- matrix(NA, (niter - nburn) / nthin, K)
    #############################################################
    
    
    ####################### Initialization ######################
    errVar <- 1 / rgamma(1, shape = a0, rate = b0)
    nu_list <- lapply(1:nmodes, function(d) {matrix(1/rgamma(pdims[d] * K, shape = 0.5, rate = 1), pdims[d], K)})
    lambda2_list <- lapply(1:nmodes, function(d) {matrix(1/rgamma(pdims[d] * K, shape = 0.5, rate = 1/as.vector(nu_list[[d]])), pdims[d], K)})
    xi_list <- 1 / rgamma(K, shape = 0.5, rate = 1)
    delta2_list <- 1 / rgamma(K, shape = 0.5, rate = 1/xi_list)
    ######### OLS/Ridge initialization ##########
    x_rows <- lapply(Xijlist, as.vector) ## Vectorize tensors X_ij's
    X_mat <- do.call(rbind, x_rows)
    gram_mat <- crossprod(X_mat)
    beta_init <- chol2inv(chol(gram_mat + diag(1e-6, nrow(gram_mat)))) %*% crossprod(X_mat, yijs)
    B <- array(beta_init, dim = pdims)
    Bfactors <- BayesCoMET:::init_CP_factors(beta_vec = beta_init, pdims = pdims, K = K)$U
    ###########################################################
    
    ###############################
    #### vectorizing the Zij arrays ####
    vecZij_list <- lapply(Zijlist,
                          function(foo) as.vector(foo)) ## vec(\mathcal{Z}_{ij})^{T}
    ###############################
    Zi_tilde_list <- lapply(1:n,
                            function(i) {
                                rows_i <- mis_starts[i]:mis_cumsum[i]
                                do.call("rbind", vecZij_list[rows_i])
                            })
    cts <- 0
    startTime <- proc.time()
    for(its in 1:niter) {
        if(its %% 1000 == 0) cat("iteration: ", its, "\n")
        cycle2Samp <- sampleB_SG(yijs = yijs, Xijlist = Xijlist,
                                 Zi_tilde_list = Zi_tilde_list,
                                 mis = mis,
                                 B_list = Bfactors, K = K,
                                 a0 = a0, b0 = b0, errVar = errVar,
                                 lambda2_list = lambda2_list, delta2_list = delta2_list,
                                 nu_list = nu_list, xi_list = xi_list,
                                 Sigma = Sigma)
        
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
        
        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            for(d in 1:nmodes) {
                lambda2Samplist[[d]][cts, ] <- as.vector(cycle2Samp$lambda2Samplist[[d]])
            }
            betaSamp[cts, ] <- as.vector(cycle2Samp$BSamp)
            errVarSamp[cts] <- cycle2Samp$errVarSamp
            delta2Samp[cts, ] <- cycle2Samp$delta2Samplist
        }
    }
    endTime <- proc.time()
    
    list(betaSamp = betaSamp, errVarSamp = errVarSamp,
         lambda2Samplist = lambda2Samplist, delta2Samp = delta2Samp,
         sampler_time = endTime - startTime)
}


post_pred_newsubj <- function(betaSamp, errVarSamp, Sigma,
                              xlist_test, zlist_test, mis, nom.level) {

    N_test <- sum(mis); n_test <- length(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)

    resids <- vector("list", n_test)
    

    #### vectorizing the Zij arrays ####
    vecZij_list_test <- lapply(dat_test$Zijlist,
                               function(foo) as.vector(foo)) ## vec(\mathcal{Z}_{ij})^{T}
    ###############################
    Zi_tilde_list_test <- lapply(1:n_test, function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        do.call("rbind", vecZij_list_test[rows_i])
    })

    ylist_test <- lapply(1:n_test, function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        dat_test$yijs[rows_i]
    })
    vecXij_test <- lapply(dat_test$Xijlist, function(foo) {as.vector(foo)})
    Xi_test <- lapply(1:n_test, function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        do.call("rbind", vecXij_test[rows_i])
    })

    for (gg in 1:n_test) {
        yhat <- list()
        for (tt in 1:length(errVarSamp)) {
            predCov <- errVarSamp[tt] * (Zi_tilde_list_test[[gg]] %*% tcrossprod(Sigma, Zi_tilde_list_test[[gg]]) + diag(1, mis[gg]))
            yhat[[tt]] <- Xi_test[[gg]] %*% betaSamp[tt, ] + drop(crossprod(chol(predCov), rnorm(mis[gg])))
        }

        yhat_samples[[gg]] <- do.call(cbind, yhat)
        
        preds[[gg]] <- rowMeans(yhat_samples[[gg]])
        #resids[[gg]] <- ylist_test[[gg]] - preds[[gg]]

        grpvec <- 1:mis[gg]

        qlower[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, (1 - nom.level) / 2))
        qupper[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, 1 - (1 - nom.level) / 2))

    }
    preds <- unlist(preds)
    mspe <- mean(unlist(lapply(resids, function(foo) {mean(foo ^ 2)})))

    covgs_95 <- unlist(covgs_95)
    ciwidth_95 <- unlist(ciwidth_95)

    covgs_90 <- unlist(covgs_90)
    ciwidth_90 <- unlist(ciwidth_90)

    covgs_80 <- unlist(covgs_80)
    ciwidth_80 <- unlist(ciwidth_80)

    res_Pred <- list(preds = preds, mspe = mspe,
                     ytrue = unlist(ylist_test),
                     covgs_95 = covgs_95, ciwidth_95_summary = summary(ciwidth_95),
                     covgs_90 = covgs_90, ciwidth_90_summary = summary(ciwidth_90),
                     covgs_80 = covgs_80, ciwidth_80_summary = summary(ciwidth_80),
                     q025 = unlist(q025), q975 = unlist(q975),
                     q05 = unlist(q05), q95 = unlist(q95),
                     q10 = unlist(q10), q90 = unlist(q90), y1hatmat = y1hatmat) #,
#                     )

    res_Pred
}


post_pred_existingsubj <- function(sampler_res, Lkron,
                                   Xijlist_test, Zijlist_test, mis, nom.level) {
    
    betaSamp <- sampler_res$betaSamp
    errVarSamp <- sampler_res$errVarSamp
    niter <- length(errVarSamp)
    N_test <- sum(mis); n_test <- length(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    
    resids <- vector("list", n_test); preds <- vector("list", n_test)
    yhat_samples <- vector("list", n_test)
    qlower <- vector("list", n_test); qupper <- vector("list", n_test)
    
    #### vectorizing the Zij arrays ####
    vecZij_list_test <- lapply(Zijlist_test,
                               function(foo) as.vector(foo)) ## vec(\mathcal{Z}_{ij})^{T}
    ###############################
    Zi_tilde_list_test <- lapply(1:n_test, function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        do.call("rbind", vecZij_list_test[rows_i])
    })
    
    vecXij_test <- lapply(Xijlist_test, function(foo) {as.vector(foo)})
    Xi_test <- lapply(1:n_test, function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        do.call("rbind", vecXij_test[rows_i])
    })
    
    for (gg in 1:n_test) {
        yhat <- list()
        for (tt in 1:niter) {
            bi <- rnorm(ncol(Zi_tilde_list_test[[gg]]), 0, sqrt(errVarSamp)) 
            ranComp <- drop(Zi_tilde_list_test[[gg]] %*% (Lkron %*% bi))
            mu_y <- drop(Xi_test[[gg]] %*% betaSamp[tt, ]) + ranComp 
            yhat[[tt]] <- rnorm(mis[gg], mean = mu_y, sd = sqrt(errVarSamp[tt]))
        }
        
        yhat_samples[[gg]] <- do.call(cbind, yhat)
        preds[[gg]] <- rowMeans(yhat_samples)
        
        grpvec <- 1:mis[gg]
        
        qlower[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, (1 - nom.level) / 2))
        qupper[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, 1 - (1 - nom.level) / 2))
    }
    
    result <- list(ypred = preds, yhat_samples = yhat_samples,
                   lower_pi = qlower, upper_pi = qupper)
    result
}

