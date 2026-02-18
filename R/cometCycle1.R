cometCycle1 <- function(yijs, xijlist, comp_zijlist, zi_tilde_list, mis,
                         B, errVar,
                         Gamma_list, Sigma_gamma_list, kdims, RRtkron) {

    n <- length(mis); N <- sum(mis)
    # kdims <- unlist(lapply(R_list, nrow))
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    nmodes <- length(kdims)

    #### calculate residuals given B
    check_y <- sapply(1:N, function(ij) {yijs[ij] - sum(as.vector(xijlist[[ij]]) * as.vector(B))})
    check_ylist <- lapply(1:n, function(i) {rows_i <- mis_starts[i]:mis_cumsum[i]; check_y[rows_i]})

    #### Impute missing data \tilde D_i's ####
    di <- vector("list", n)
    Di_tilde <- vector("list", n)
    Gkron <- revkronAll(Gamma_list)
    for (i in 1:n) {
        Vyd <- errVar * zi_tilde_list[[i]] %*% Gkron %*% RRtkron
        VyyInv <- chol2inv(chol(Vyd %*% crossprod(Gkron, t(zi_tilde_list[[i]])) + errVar * diag(1, nrow(zi_tilde_list[[i]]))))
        VydTvyyInv <- crossprod(Vyd, VyyInv)
        postRanVar <- errVar * RRtkron - drop(VydTvyyInv %*% Vyd)
        postRanMean <- drop(VydTvyyInv %*% check_ylist[[i]])
        di[[i]] <- postRanMean + drop(crossprod(chol(postRanVar + diag(1e-6, nrow(postRanVar))), rnorm(length(postRanMean))))
        Di_tilde[[i]] <- array(di[[i]], dim = kdims)
    }
    #message("Imputation of di's Done")
    for(d in 1:nmodes) {
        check_z_gamma_groups <- vector("list", n)
        for(i in 1:n) {
            rows_i <- mis_starts[i]:mis_cumsum[i]
            check_z_gamma_groups[[i]]  <- do.call("rbind", lapply(comp_zijlist[rows_i],
                                                                  function(foo) {as.vector(mode_matrix(foo, d) %*% tcrossprod(LOOrevkron(Gamma_list, d), mode_matrix(Di_tilde[[i]], d)))}))
        }
        check_z_gamma <- do.call("rbind", check_z_gamma_groups)

        postVar_gamma <- chol2inv(chol(crossprod(check_z_gamma) / errVar +
                                           chol2inv(chol(Sigma_gamma_list[[d]]))))
        postMean_gamma <- drop(postVar_gamma %*% crossprod(check_z_gamma, check_y) / errVar)

        gammaSamp <- postMean_gamma + drop(crossprod(chol(postVar_gamma + diag(1e-6, nrow(postVar_gamma))), rnorm(length(postMean_gamma))))
        #message("Done Gamma", d)
        Gamma_list[[d]] <- matrix(gammaSamp, kdims[d], kdims[d])
    }

    list(GammaSamplist = Gamma_list, Di_tilde = Di_tilde)
}
