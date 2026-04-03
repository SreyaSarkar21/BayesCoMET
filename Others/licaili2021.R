library(Matrix)
library(glmnet)
library(scalreg)
library(parallel)

CV_fixed_effects <- function(xlist, ylist, zlist, mi_train, a.seq=seq(0.5, 10, 0.5), inf.coord,
                              standardize, intercept, std_response){
    n <- length(mi_train)
    n.tr <- round(0.75*n) # 75% subjects to train and rest to test
    best.err <- sum(unlist(ylist)^2)
    best.a <- a.seq[1]

    for(ia in 1:length(a.seq)){
        est.re <- Fixed_effects_estimation2(xlist[1:n.tr], ylist[1:n.tr],
                                            zlist[1:n.tr], mi_train = mi_train[1:n.tr],
                                            a = a.seq[ia], inf.coord = inf.coord,
                                            standardize = standardize, intercept = intercept,
                                            std_response = std_response)

        pred.err <- sum((unlist(ylist[-(1:n.tr)]) - drop(do.call(rbind, xlist[-(1:n.tr)]) %*% est.re$beta.hat.lasso))^2) # For intercept = FALSE
        if(pred.err < best.err){
            best.err <- pred.err
            best.a <- a.seq[ia]
        }
    }
    list(best.a = best.a)
}

Fixed_effects_estimation <- function(xlist, ylist, zlist, mi_train, a, inf.coord = NULL,
                                      standardize, intercept, std_response){

    N <- sum(mi_train)
    n <- length(mi_train)
    mi_ends <- cumsum(mi_train)
    mi_starts <- c(1, mi_ends[-length(mi_train)] + 1)

    x <- do.call(rbind, xlist)
    y <- unlist(ylist)
    p <- ncol(x)

    ## mixed effect model fitting
    x.a <- x
    y.a <- y
    tr.a <- 0

    for (i in 1:n){
        cur.mem <- mi_starts[i]:mi_ends[i]
        sigmai <- a * tcrossprod(zlist[[i]]) + diag(nrow(zlist[[i]]))
        Ri <- chol(sigmai)
        Sig.a.inv.half <- backsolve(Ri, diag(nrow(Ri)))
        tr.a <- tr.a + sum(Sig.a.inv.half^2)

        x.a[cur.mem, ] <- Sig.a.inv.half %*% x[cur.mem, ]
        y.a[cur.mem] <- Sig.a.inv.half %*% as.matrix(y[cur.mem])
    }

    sig.init <- scalreg::scalreg(x.a, y.a)$hsigma

    lasso.fit <- glmnet::glmnet(
        x.a, y.a,
        lambda = sig.init * sqrt(2 * log(p) / N),
        standardize = standardize,
        intercept = intercept,
        standardize.response = std_response
    )

    beta.hat.lasso <- as.vector(glmnet::coef.glmnet(lasso.fit))[-1]

    ## initialize full-length vectors
    beta.hat.db <- rep(NA_real_, p)
    beta.hat.db.sd <- rep(NA_real_, p)

    if (is.null(inf.coord)) {
        return(list(beta.hat.lasso = beta.hat.lasso,
                    beta.hat.db = beta.hat.db,
                    tr.a = tr.a))
    }

    res <- y.a - x.a %*% beta.hat.lasso

    ## use cores allocated by scheduler
    cores <- as.integer(Sys.getenv("NSLOTS", 1))
    print(cores)
    cores <- max(2, cores)

    res_db_list <- parallel::mclapply(
        inf.coord,
        function(col.j) {
            sig.x <- scalreg::scalreg(x.a[, -col.j, drop = FALSE], x.a[, col.j])$hsigma

            kappa.hat.mlm <- glmnet::glmnet(
                x.a[, -col.j, drop = FALSE], x.a[, col.j],
                lambda = sig.x * sqrt(2 * log(p) / N),
                standardize = standardize,
                intercept = intercept,
                standardize.response = std_response
            )

            gam.j <- as.vector(glmnet::coef.glmnet(kappa.hat.mlm))[-1]
            wj.mlm <- x.a[, col.j] - x.a[, -col.j, drop = FALSE] %*% gam.j

            denom <- drop(crossprod(wj.mlm, x.a[, col.j]))
            if (!is.finite(denom) || abs(denom) < 1e-12) {
                return(c(beta = NA_real_, sd = NA_real_))
            }

            beta.hat.db.j <- beta.hat.lasso[col.j] +
                drop(crossprod(wj.mlm, res)) / denom

            num <- 0
            for (i in 1:n) {
                cur.mem <- mi_starts[i]:mi_ends[i]
                num <- num + (sum(wj.mlm[cur.mem] * res[cur.mem]))^2
            }

            beta.hat.db.sd.j <- sqrt(num) / denom
            c(beta = beta.hat.db.j, sd = beta.hat.db.sd.j)
        },
        mc.cores = cores,
        mc.set.seed = TRUE
    )

    res_db_mat <- do.call(rbind, res_db_list)
    colnames(res_db_mat) <- c("beta", "sd")

    beta.hat.db[inf.coord] <- res_db_mat[, "beta"]
    beta.hat.db.sd[inf.coord] <- res_db_mat[, "sd"]

    list(beta.hat.lasso = beta.hat.lasso,
         beta.hat.db = beta.hat.db,
         beta.hat.db.sd = beta.hat.db.sd,
         tr.a = tr.a,
         xmean = colMeans(x),
         xsd = apply(x, 2, sd))
}

##### For conditional prediction
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
