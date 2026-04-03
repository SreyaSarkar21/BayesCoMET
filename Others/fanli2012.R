############### Fan and Li (2012) ###############
######## For Simulation Study ########
######## Function to calculate the number of times the correct fixed effects are selected ########


calc_CF <- function(ylist, xlist, zlist,
                    capM, beta.true, xtrain_sds,
                    standardize, intercept, std_response,
                    newylist, newxlist, newzlist) {
    library(Matrix)
    library(glmnet)
    zmat <- as.matrix(bdiag(zlist)) # N x nq
    xmat <- do.call(rbind, xlist) # N x p
    y <- matrix(unlist(ylist), ncol = 1) # N x 1

    Pz <- chol2inv(chol(diag(length(y)) + tcrossprod(zmat %*% capM, zmat)))
    Pz_svd <- svd(Pz)
    Pz_half <- Pz_svd$u %*% tcrossprod(diag(sqrt(Pz_svd$d)), Pz_svd$u)
    y_scaled <- Pz_half %*% y
    xmat_scaled <- Pz_half %*% xmat

    cv.out <- cv.glmnet(xmat_scaled, y_scaled, alpha = 1,
                        intercept = intercept, standardize = standardize, standardize.response = std_response)

    # lambda value chosen by 10-fold CV
    beta.hat.lmin <- coef(cv.out, s = cv.out$lambda.min)[-1]

    rmse.lmin <- sqrt(mean((beta.hat.lmin - beta.true) ^ 2))
    beta.hat.lmin.orig.scale <- beta.hat.lmin / xtrain_sds
    rmse.lmin.scale <- sqrt(mean((beta.hat.lmin.orig.scale - beta.true) ^ 2))

    true_indx.lmin <- which(beta.true != 0)
    est_indx.lmin <- which(beta.hat.lmin != 0)
    cf.lmin <- length(intersect(est_indx.lmin, true_indx.lmin))
    tpr.lmin <- cf.lmin / length(true_indx.lmin)
    fpr.lmin <- length(intersect(est_indx.lmin, which(beta.true == 0))) / length(which(beta.true == 0))

    ## Prediction (Marginal)

    y.test <- as.matrix(unlist(newylist), ncol = 1)
    x.test <- do.call(rbind, newxlist)
    y.hat.lmin <- x.test %*% beta.hat.lmin

    mspe.lmin <- mean((y.hat.lmin - y.test) ^ 2)

    list(beta.hat.lmin = beta.hat.lmin, cf.lmin = cf.lmin,
         tpr.lmin = tpr.lmin, fpr.lmin = fpr.lmin,
         rmse.lmin = rmse.lmin,
         y.hat.lmin = y.hat.lmin, ytrue.test = y.test, mspe.lmin = mspe.lmin,
         beta.hat.lmin.orig.scale = beta.hat.lmin.orig.scale,
         rmse.lmin.scale = rmse.lmin.scale,
         xmat_scaled = xmat_scaled, y_scaled = y_scaled,
         Pz = Pz, Pz_half = Pz_half)
}


######## Function to calculate coverage of confidence intervals for beta ########
########(following Debiased estimator by LiCaiLi2021) ########

covg_fanli <- function(lasso_est, xtrain_sds, beta.true, n, m, intercept, standardize, std_response) {
    library(Matrix)
    library(glmnet)
    p <- length(beta.true)
    y_scaled <- lasso_est$y_scaled
    xmat_scaled <- lasso_est$xmat_scaled

    beta.hat.db.sd.lmin <- rep(NA, p)
    beta.hat.db.lmin <- rep(NA, p)

    res.lmin <- y_scaled - xmat_scaled %*% lasso_est$beta.hat.lmin

    for(j in 1:p){
        kappa.hat.mlm <- cv.glmnet(as.matrix(xmat_scaled[, -j], ncol = p-1), as.vector(xmat_scaled[, j]),
                                   alpha = 1, nfolds = 10,
                                   intercept = intercept, standardize = standardize, standardize.response = std_response)
        gam.j.lmin <- coef(kappa.hat.mlm, s = kappa.hat.mlm$lambda.min)[-1]
        wj.mlm.lmin <- xmat_scaled[, j] - xmat_scaled[, -j] %*% gam.j.lmin
        beta.hat.db.lmin[j] <- lasso_est$beta.hat.lmin[j] + sum(wj.mlm.lmin * res.lmin)/sum(wj.mlm.lmin * xmat_scaled[,j])

        num.lmin <- 0
        for(i in 1:n){
            cur.mem <- ((i-1)*m+1):(i*m)
            num.lmin <- num.lmin + sum(wj.mlm.lmin[cur.mem] * res.lmin[cur.mem])^2
        }

        beta.hat.db.sd.lmin[j] <- sqrt(num.lmin)/sum(wj.mlm.lmin * as.vector(xmat_scaled[,j]))
    }

    beta.hat.db.lmin.scaled <- beta.hat.db.lmin / xtrain_sds
    beta.hat.db.sd.lmin.scaled <- beta.hat.db.sd.lmin / xtrain_sds
    z_95 <- qnorm(0.975); z_90 <- qnorm(0.95); z_80 <- qnorm(0.90)

    covgs_95.lmin <- rep(NA, p)
    covgs_90.lmin <- rep(NA, p)
    covgs_80.lmin <- rep(NA, p)

    for(pp in 1:p) {
        covgs_95.lmin[pp] <- as.numeric(beta.hat.db.lmin[pp] - z_95 * beta.hat.db.sd.lmin[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin[pp] + z_95 * beta.hat.db.sd.lmin[pp])
        covgs_90.lmin[pp] <- as.numeric(beta.hat.db.lmin[pp] - z_90 * beta.hat.db.sd.lmin[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin[pp] + z_90 * beta.hat.db.sd.lmin[pp])
        covgs_80.lmin[pp] <- as.numeric(beta.hat.db.lmin[pp] - z_80 * beta.hat.db.sd.lmin[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin[pp] + z_80 * beta.hat.db.sd.lmin[pp])
    }

    covgs_95.lmin.scaled <- rep(NA, p)
    covgs_90.lmin.scaled <- rep(NA, p)
    covgs_80.lmin.scaled <- rep(NA, p)

    for(pp in 1:p) {
        covgs_95.lmin.scaled[pp] <- as.numeric(beta.hat.db.lmin.scaled[pp] - z_95 * beta.hat.db.sd.lmin.scaled[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin.scaled[pp] + z_95 * beta.hat.db.sd.lmin.scaled[pp])
        covgs_90.lmin.scaled[pp] <- as.numeric(beta.hat.db.lmin.scaled[pp] - z_90 * beta.hat.db.sd.lmin.scaled[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin.scaled[pp] + z_90 * beta.hat.db.sd.lmin.scaled[pp])
        covgs_80.lmin.scaled[pp] <- as.numeric(beta.hat.db.lmin.scaled[pp] - z_80 * beta.hat.db.sd.lmin.scaled[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin.scaled[pp] + z_80 * beta.hat.db.sd.lmin.scaled[pp])
    }

    list(beta.hat.db.lmin = beta.hat.db.lmin,
         beta.hat.db.sd.lmin = beta.hat.db.sd.lmin,
         covgs_95.lmin = covgs_95.lmin, covgs_90.lmin = covgs_90.lmin, covgs_80.lmin = covgs_80.lmin,
         beta.hat.db.lmin.scaled = beta.hat.db.lmin.scaled,
         beta.hat.db.sd.lmin.scaled = beta.hat.db.sd.lmin.scaled,
         covgs_95.lmin.scaled = covgs_95.lmin.scaled, covgs_90.lmin.scaled = covgs_90.lmin.scaled, covgs_80.lmin.scaled = covgs_80.lmin.scaled)
}


#### For Conditional Prediction

compute_blup_quasi <- function(resids, Z, sigma_Z, grp) {
    grp_uniq <- unique(grp)

    n_dim <- length(grp_uniq)
    m_dim <- length(grp)
    stopifnot(nrow(Z) == m_dim)

    q_dim <- ncol(Z)
    blup  <- matrix(NA_real_, nrow = q_dim, ncol = n_dim)

    for (gx in seq_len(n_dim)) {
        ixs <- which(grp == grp_uniq[gx])
        r_g <- resids[ixs]
        Z_g <- Z[ixs, , drop = FALSE]

        # Sigma_Y_g = Z_g Sigma_Z Z_g^T + I
        SigmaY_g <- tcrossprod(Z_g %*% sigma_Z, Z_g) + diag(1, nrow(Z_g))

        # b_g = Sigma_Z Z_g^T SigmaY_g^{-1} r_g
        blup[, gx] <- sigma_Z %*% crossprod(Z_g, solve(SigmaY_g, r_g))
    }

    list(group_ids = grp_uniq, blup = blup)
}


predict_quasi <- function(X_test, Z_test, b_hat, blup_map, grp_test) {
    stopifnot(nrow(X_test) == nrow(Z_test))
    m_dim <- nrow(X_test)

    # baseline (marginal / population-level) prediction
    y_hat <- as.numeric(X_test %*% b_hat)

    # map each test row's group to a BLUP column (NA if unseen group)
    idx <- match(grp_test, blup_map$group_ids)

    has <- which(!is.na(idx))
    if (length(has) > 0) {
        # add Z * b_i rowwise
        y_hat[has] <- y_hat[has] + rowSums(Z_test[has, , drop = FALSE] *
                                               t(blup_map$blup[, idx[has], drop = FALSE]))
    }

    y_hat
}

########### For Data Analysis ############
fanli_pred_data <- function(ylist, xlist, zlist, standardize, intercept, std_response,
                          newylist, newxlist, newzlist) {
    library(Matrix)
    library(glmnet)
    library(MASS)

    mi_train <- unlist(lapply(ylist, length)); n_train <- length(ylist)
    mi_test  <- unlist(lapply(newylist, length)); n_test <- length(newylist)

    xmat    <- do.call(rbind, xlist) # N × p
    y_train <- matrix(unlist(ylist), ncol=1) # N × 1


    ## Extracting scalar from capM for memory-safe block-wise computation of Pz
    cval <- log(sum(mi_train))

    ## Pz_list: list of per-subject (I_{mi} + c Zi Zi^T)^{-1}
    Pz_list <- vector("list", length(zlist))

    for (i in seq_along(zlist)) {
        Zi <- as.matrix(zlist[[i]]) # mi × q
        mi <- nrow(Zi)

        ## small mi × mi matrix: A_i = I + c * Zi Zi^T
        A_i <- diag(mi) + cval * (Zi %*% t(Zi))

        ## inverse via Cholesky
        Pz_list[[i]] <- chol2inv(chol(A_i))
    }

    ## assemble full Pz as block diagonal WITHOUT dense expansion
    Pz <- as.matrix(bdiag(Pz_list))

    ## For each block, compute symmetric square-root of Pz_i
    Pz_half_list <- lapply(Pz_list, function(Pi) {
        eig <- eigen(Pi, symmetric = TRUE)
        eig$vectors %*% (sqrt(eig$values) * t(eig$vectors))
    })

    ## block diagonal concatenation
    Pz_half <- as.matrix(bdiag(Pz_half_list))

    ## Scaling y and x using Pz_half
    y_scaled <- Pz_half %*% y_train
    xmat_scaled <- Pz_half %*% xmat

    cv.result <- cv.glmnet(xmat_scaled, y_scaled, alpha = 1,
                           standardize = standardize,
                           intercept = intercept,
                           standardize.response = std_response)

    beta.hat.lmin <- as.vector(coef(cv.result, s = cv.result$lambda.min))

    ## predictions on TRAINING (scaled y)
    pred_y_scaled_lmin <- cbind(rep(1, nrow(xmat_scaled)), xmat_scaled) %*% beta.hat.lmin

    ## Out-of-Sample Prediction (Marginal)
    y.test <- as.matrix(unlist(newylist), ncol = 1)
    x.test <- cbind(rep(1, nrow(y.test)), do.call(rbind, newxlist))

    y.hat.lmin <- drop(x.test %*% beta.hat.lmin)

    ## group-wise weighted MSPE and MAD
    groups <- rep(seq_along(mi_test), times = mi_test)
    ## lambda_min
    resids_list <- split(y.hat.lmin - drop(y.test), groups)
    mspe.lmin <- mean(sapply(resids_list, function(foo) mean(foo^2)))

    list(beta.hat.lmin = beta.hat.lmin,
         xmat_scaled = xmat_scaled, y_scaled = y_scaled,
         Pz = Pz, Pz_half = Pz_half,
         cv.result = cv.result,
         pred_y_scaled_lmin = pred_y_scaled_lmin,
         y.hat.lmin = y.hat.lmin,
         y_true = drop(y.test), x.test = x.test,
         mspe.lmin = mspe.lmin, cval = cval)
}


