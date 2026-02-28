

generate_cpB_with_sparse_factors <- function(pdims, K = 2,
                                                   sparsity = 0.25,
                                                   seed = 123) {
    set.seed(seed)
    nmodes <- length(pdims)

    if(nmodes == 2) {
        stopifnot(K <= min(pdims))
    }

    # Generate sparse factor matrices
    factor_list <- vector("list", nmodes)
    for (d in 1:nmodes) {
        total_entries <- pdims[d] * K
        vec <- rep(0, total_entries)
        n_nonzero <- round(sparsity * total_entries)
        vec[sample(total_entries, n_nonzero)] <- sample(c(-2, -1, 1, 2), n_nonzero, replace = TRUE)
        factor_list[[d]] <- Matrix::Matrix(vec, nrow = pdims[d], ncol = K)
    }

    # Construct CP decomposition: sum of K rank-1 tensors
    if(nmodes == 2) {
        B <- tcrossprod(factor_list[[1]], factor_list[[2]])
    } else {
        # Initialize the tensor
        B <- array(0, dim = pdims)

        # Construct CP decomposition: sum of K rank-1 tensors
        for (k in 1:K) {
            # Start with first vector
            rank1 <- factor_list[[1]][, k]
            for (d in 2:nmodes) {
                rank1 <- outer(rank1, factor_list[[d]][, k])
                dim(rank1) <- c(dim(rank1))  # reshape for next step
            }
            B <- B + array(rank1, dim = pdims)
        }
    }

    return(list(B = B, B_factors = factor_list))
}

genarrayNormal <- function(dims, Mean, covlist) {
    covWhole = revkronAll(lapply(covlist, chol))
    vecX <- as.vector(Mean) + drop(rnorm(prod(dims)) %*% covWhole)
    X <- array(vecX, dim = dims)
}


simdata <- function(pdims, qdims, n, m, errVar, B, L_list,
                    xcov_var, zcov_var, myseed) {
    set.seed(myseed)

    nmodes <- length(pdims)
    sdims <- unlist(lapply(L_list, ncol))

    N <- n * m
    group <- rep(1:n, each = m)

    #### generate covariates from N(0, xcov_var * I) and N(0, zcov_var * I) ####
    Xijlist <- vector("list", N); Zijlist <- vector("list", N)
    for(ij in 1:N) {
        Xijlist[[ij]] <- genarrayNormal(dims = pdims,
                                        Mean = array(0, dim = pdims),
                                        covlist = lapply(pdims, function(foo) {xcov_var * diag(foo)}))
        Zijlist[[ij]] <- genarrayNormal(dims = qdims,
                                        Mean = array(0, dim = qdims),
                                        covlist = lapply(qdims, function(foo) {zcov_var * diag(foo)}))
    }

    #### generate response ####
    yijs <- numeric(length = N)
    for(i in 1:n) {
        vecDi <- matrix(rnorm(prod(sdims)), ncol = 1)
        vecAi <- drop(revkronAll(L_list) %*% vecDi)
        for(j in ((i-1)*m+1):(i*m)) {
            fixef_comp <- sum(as.vector(Xijlist[[j]]) * as.vector(B))
            ranef_comp <- sum(vecAi * as.vector(Zijlist[[j]]))
            yijs[j] <- fixef_comp + ranef_comp + rnorm(1, 0, sd = sqrt(errVar))
        }
    }

    list(yijs = yijs, Xijlist = Xijlist, Zijlist = Zijlist,
         B = B, L_list = L_list, errVar = errVar)
}
