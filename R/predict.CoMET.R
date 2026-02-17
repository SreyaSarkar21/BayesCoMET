#' @title predict.CoMET
#'
#' @description This function performs posterior predictive sampling using the compressed mixed-effects tensor (CoMET) model.
#'
#' @param sampler_res a list of posterior samples of \eqn{(\beta, \tau^2, \gamma)} returned by the function \code{\link{CoMET}}.
#' @param ylist_test a list of response vectors where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param xlist_test a list of fixed effect covariates where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param zlist_test a list of random effect covariates where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param R_list the \eqn{k_1 \times q} random projection matrices used to implement the collapsed Gibbs sampler using \code{\link{CoMET}}.
#' @param S_list the \eqn{k_2 \times q} random projection matrices used to implement the collapsed Gibbs sampler using \code{\link{CoMET}}.
#' @param nom.level nominal level for constructing prediction intervals.
#' @return a list containing the following components:
#' \describe{
#' \item{ypred}{predicted values of response.}
#' \item{mspe}{mean square prediction error.}
#' \item{lower_pi}{lower limits of the prediction intervals.}
#' \item{upper_pi}{upper limits of the prediction intervals.}
#' \item{yhat_samples}{posterior samples of the response.}
#' }
#' @importFrom stats rnorm rgamma quantile
#' @export

predict.CoMET <- function(betaSamp, errVarSamp, gammaSamplist, R_list, S_list,
                          yijs_test, Xijlist_test, Zijlist_test, mis,
                          burnin, thin) {
    
    nmodes <- length(gammaSamplist)
    betaSamp <- betaSamp[seq(burnin + 1, nrow(betaSamp), by = thin), ]
    errVarSamp <- errVarSamp[seq(burnin + 1, length(errVarSamp), by = thin)]
    gammaSamplist <- lapply(seq_along(gammaSamplist),
                            function(d) {gammaSamplist[[d]][seq(burnin + 1, nrow(gammaSamplist[[d]]), by = thin), ]})
    kdims <- unlist(lapply(R_list, nrow))
    n_test <- length(mis); N_test <- sum(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    
    resids_mean <- vector("list", n_test)
    resids_median <- vector("list", n_test)
    preds_mean <- vector("list", n_test)
    preds_median <- vector("list", n_test)
    q025 <- vector("list", n_test); q975 <- vector("list", n_test)
    covgs_95 <- vector("list", n_test); ciwidth_95 <- vector("list", n_test)
    q05 <- vector("list", n_test); q95 <- vector("list", n_test)
    covgs_90 <- vector("list", n_test); ciwidth_90 <- vector("list", n_test)
    q10 <- vector("list", n_test); q90 <- vector("list", n_test)
    covgs_80 <- vector("list", n_test); ciwidth_80 <- vector("list", n_test)
    y1hatmat <- matrix(NA, mis[1], length(errVarSamp))
    
    
    #### compress Zij arrays ####
    comp_Zijlist_test <- vector("list", N_test)
    Skron_not1 <- LOORevKron(S_list, d = 1)
    for(ij in 1:N_test) {
        Zij_tilde_mode1 <- S_list[[1]] %*% tcrossprod(mode_matrix(Zijlist_test[[ij]], n = 1), Skron_not1)
        comp_Zijlist_test[[ij]] <- array(Zij_tilde_mode1, dim = kdims)
    }
    #### vectorize the compressed Zij arrays ####
    vec_comp_Zijlist_test <- lapply(comp_Zijlist_test,
                                    function(foo) {as.vector(foo)})
    ###############################
    Zi_tilde_list_test <- lapply(1:n_test,  function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        do.call("rbind", vec_comp_Zijlist_test[rows_i])
    })
    
    ylist_test <- lapply(1:n_test,
                         function(i) {
                             rows_i <- mis_starts[i]:mis_cumsum[i]
                             yijs_test[rows_i]})
    vecXij_test <- lapply(Xijlist_test, function(foo) {as.vector(foo)})
    Xi_test <- lapply(1:n_test,
                      function(i) {
                          rows_i <- mis_starts[i]:mis_cumsum[i]
                          do.call("rbind", vecXij_test[rows_i])
                      })
    
    GammaSamplist <- lapply(seq_len(nmodes),
                            function(d) {lapply(seq_len(nrow(gammaSamplist[[d]])),
                                                function(foo) matrix(gammaSamplist[[d]][foo, ], kdims[d], kdims[d]))})
    
    for (gg in 1:n_test) {
        yhat <- list()
        for (tt in 1:length(errVarSamp)) {
            obj <- lapply(seq_len(nmodes), function(d) GammaSamplist[[d]][[tt]])
            Gkron <- kronAll(obj)
            Rkron <- kronAll(lapply(R_list, tcrossprod))
            predCov <- errVarSamp[tt] * (Zi_tilde_list_test[[gg]] %*% Gkron %*% Rkron %*% crossprod(Gkron, t(Zi_tilde_list_test[[gg]])) + diag(1, mis[gg]))
            yhat[[tt]] <- Xi_test[[gg]] %*% betaSamp[tt, ] + drop(crossprod(chol(predCov), rnorm(mis[gg])))
        }
        
        yhatmat <- do.call(cbind, yhat)
        if(gg == 1) {
            y1hatmat <- yhatmat
        }
        
        preds_mean[[gg]] <- rowMeans(yhatmat)
        resids_mean[[gg]] <- ylist_test[[gg]] - preds_mean[[gg]]
        
        preds_median[[gg]] <- apply(yhatmat, 1, median)
        resids_median[[gg]] <- ylist_test[[gg]] - preds_median[[gg]]
        
        grpvec <- 1:length(ylist_test[[gg]])
        
        q025[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.025))
        q975[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.975))
        covgs_95[[gg]] <- sapply(grpvec, function(ii) {as.numeric(q025[[gg]][ii] <= ylist_test[[gg]][ii] & ylist_test[[gg]][ii] <= q975[[gg]][ii])})
        ciwidth_95[[gg]] <- sapply(grpvec, function(ii) q975[[gg]][ii] - q025[[gg]][ii])
        
        q05[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.05))
        q95[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.95))
        covgs_90[[gg]] <- sapply(grpvec, function(ii) {as.numeric(q05[[gg]][ii] <= ylist_test[[gg]][ii] & ylist_test[[gg]][ii] <= q95[[gg]][ii])})
        ciwidth_90[[gg]] <- sapply(grpvec, function(ii) q95[[gg]][ii] - q05[[gg]][ii])
        
        q10[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.10))
        q90[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.90))
        covgs_80[[gg]] <- sapply(grpvec, function(ii) {as.numeric(q10[[gg]][ii] <= ylist_test[[gg]][ii] & ylist_test[[gg]][ii] <= q90[[gg]][ii])})
        ciwidth_80[[gg]] <- sapply(grpvec, function(ii) q90[[gg]][ii] - q10[[gg]][ii])
        
    }
    #preds_mean <- unlist(preds_mean); preds_median <- unlist(preds_median)
    mspe_mean <- mean(unlist(lapply(resids_mean, function(foo) {mean(foo ^ 2)})))
    mspe_median <- mean(unlist(lapply(resids_median, function(foo) {mean(foo ^ 2)})))
    mad_mean <- mean(unlist(lapply(resids_mean, function(foo) {mean(abs(foo))})))
    mad_median <- mean(unlist(lapply(resids_median, function(foo) {mean(abs(foo))})))
    
    covgs_95 <- unlist(covgs_95)
    ciwidth_95 <- unlist(ciwidth_95)
    
    covgs_90 <- unlist(covgs_90)
    ciwidth_90 <- unlist(ciwidth_90)
    
    covgs_80 <- unlist(covgs_80)
    ciwidth_80 <- unlist(ciwidth_80)
    
    result <- list(preds_mean = preds_mean, mspe_mean = mspe_mean, mad_mean = mad_mean,
                   preds_median = preds_mean, mspe_median = mspe_median, mad_median = mad_median,
                   ytrue = yijs_test,
                   covgs_95 = covgs_95, ciwidth_95_summary = summary(ciwidth_95),
                   covgs_90 = covgs_90, ciwidth_90_summary = summary(ciwidth_90),
                   covgs_80 = covgs_80, ciwidth_80_summary = summary(ciwidth_80),
                   q025 = unlist(q025), q975 = unlist(q975),
                   q05 = unlist(q05), q95 = unlist(q95),
                   q10 = unlist(q10), q90 = unlist(q90), mi_test = mis,
                   y1hatmat = y1hatmat, R_list = R_list, S_list = S_list)
    #y1hatmat = y1hatmat, R_list = R_list, S_list = S_list
    result
}

