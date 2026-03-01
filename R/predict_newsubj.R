#' @title predict_newsubj
#'
#' @description This function performs posterior predictive sampling using the CoMET model.
#'
#' @param object a list of posterior samples of parameters returned by the function \code{\link{comet}}.
#' @param kdims a vector of length \eqn{D} (number of tensor modes), where each element equals the \eqn{d}-th mode-specific compressed covariance dimension \eqn{k_d}.
#' @param R_list the list of \eqn{k_d \times q} random projection matrices used for compressing the random slope in the Gibbs sampler of CoMET. The list should be same as that used as argument to the function \code{\link{comet}}.
#' @param S_list the list of \eqn{k_d \times q} random projection matrices used for compressing the random-effect covariates in the Gibbs sampler of CoMET. The list should be same as that used as argument to the function \code{\link{comet}}.
#' @param y_test a list of response vectors where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param xlist_test a list of fixed effect covariates where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param zlist_test a list of random effect covariates where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param mis a vector of cluster sizes of the new subjects.
#' @param nom.level nominal level for constructing prediction intervals.
#' @return a list containing the following components:
#' \describe{
#' \item{ypred}{predicted values of response.}
#' \item{yhat_samples}{posterior samples of the response.}
#' \item{lower_pi}{lower limits of the prediction intervals.}
#' \item{upper_pi}{upper limits of the prediction intervals.}
#' }
#' @importFrom stats rnorm quantile
#' @export

predict_newsubj <- function(object, kdims, R_list, S_list,
                              y_test, xlist_test, zlist_test, mis, nom.level) {

    betaSamp <- object$betaSamp
    errVarSamp <- object$errVarSamp
    gammaSamplist <- object$gammaSamplist

    nmodes <- length(gammaSamplist)
    n_test <- length(mis); N_test <- sum(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    niter <- length(errVarSamp)

    preds <- vector("list", n_test)
    yhat_samples <- vector("list", n_test)
    qlower <- vector("list", n_test); qupper <- vector("list", n_test)

    #### compress Zij arrays ####
    comp_zlist_test <- vector("list", N_test)
    Skron_not1 <- revkronLOO(S_list, d = 1)
    for(ij in 1:N_test) {
        z_tilde_mode1 <- S_list[[1]] %*% tcrossprod(mode_matricize(zlist_test[[ij]], d = 1), Skron_not1)
        comp_zlist_test[[ij]] <- array(z_tilde_mode1, dim = kdims)
    }
    #### vectorize the compressed Zij arrays ####
    vec_comp_zlist_test <- lapply(comp_zlist_test,
                                    function(foo) {as.vector(foo)})
    ###############################
    z_tilde_list_test <- lapply(1:n_test,  function(i) {
        rows_i <- mis_starts[i]:mis_cumsum[i]
        do.call("rbind", vec_comp_zlist_test[rows_i])
    })

    RRtkron <- revkronAll(lapply(R_list, tcrossprod))

    vecx_test <- lapply(xlist_test, function(foo) {as.vector(foo)})

    Xi_test <- lapply(1:n_test,
                      function(i) {
                          rows_i <- mis_starts[i]:mis_cumsum[i]
                          do.call("rbind", vecx_test[rows_i])
                      })

    GammaSamplist <- lapply(seq_len(nmodes),
                            function(d) {lapply(seq_len(nrow(gammaSamplist[[d]])),
                                                function(foo) matrix(gammaSamplist[[d]][foo, ], kdims[d], kdims[d]))})

    for (gg in 1:n_test) {
        yhat <- list()
        for (tt in 1:niter) {
            obj <- lapply(seq_len(nmodes), function(d) GammaSamplist[[d]][[tt]])
            Gkron <- revkronAll(obj)
            predCov <- errVarSamp[tt] * (z_tilde_list_test[[gg]] %*% Gkron %*% RRtkron %*% crossprod(Gkron, t(z_tilde_list_test[[gg]])) + diag(1, mis[gg]))
            yhat[[tt]] <- Xi_test[[gg]] %*% betaSamp[tt, ] + drop(crossprod(chol(predCov), rnorm(mis[gg])))
        }

        yhat_samples[[gg]] <- do.call(cbind, yhat)
        preds[[gg]] <- rowMeans(yhat_samples[[gg]])

        grpvec <- 1:mis[gg]

        qlower[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, (1 - nom.level) / 2))
        qupper[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, 1 - (1 - nom.level) / 2))

    }

    result <- list(ypred = preds, yhat_samples = yhat_samples,
                   lower_pi = qlower, upper_pi = qupper)
    result
}


# predict_existingsubj <- function(object,
#                                    R_list, S_list, kdims,
#                                    y_test, xlist_test, zlist_test, mis) {
#
#     betaSamp <- object$betaSamp
#     errVarSamp <- object$errVarSamp
#     gammaSamplist <- object$gammaSamplist
#     ranefSamplist <- object$ranefSamplist
#     nmodes <- length(gammaSamplist)
#     niter <- length(errVarSamp)
#
#     n_test <- length(mis); N_test <- sum(mis)
#     mis_cumsum <- cumsum(mis)
#     mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
#
#     resids <- vector("list", n_test); preds <- vector("list", n_test)
#     yhat_samples <- vector("list", n_test)
#     qlower <- vector("list", n_test); qupper <- vector("list", n_test)
#
#
#     #### compress Zij arrays ####
#     comp_zlist_test <- vector("list", N_test)
#     Skron_not1 <- revkronLOO(S_list, d = 1)
#     for(ij in 1:N_test) {
#         z_tilde_mode1 <- S_list[[1]] %*% tcrossprod(mode_matricize(zlist_test[[ij]], d = 1), Skron_not1)
#         comp_zlist_test[[ij]] <- array(z_tilde_mode1, dim = kdims)
#     }
#     #### vectorize the compressed Zij arrays ####
#     vec_comp_zlist_test <- lapply(comp_zlist_test,
#                                     function(foo) {as.vector(foo)})
#     ###############################
#     z_tilde_list_test <- lapply(1:n_test,  function(i) {
#         rows_i <- mis_starts[i]:mis_cumsum[i]
#         do.call("rbind", vec_comp_zlist_test[rows_i])
#     })
#     # RRtkron <- BayesCoMET::revkronAll(lapply(R_list, tcrossprod))
#     ylist_test <- lapply(1:n_test,
#                          function(i) {
#                              rows_i <- mis_starts[i]:mis_cumsum[i]
#                              y_test[rows_i]})
#     vecx_test <- lapply(xlist_test, function(foo) {as.vector(foo)})
#     Xi_test <- lapply(1:n_test,
#                       function(i) {
#                           rows_i <- mis_starts[i]:mis_cumsum[i]
#                           do.call("rbind", vecx_test[rows_i])
#                       })
#
#     GammaSamplist <- lapply(seq_len(nmodes),
#                             function(d) {lapply(seq_len(nrow(gammaSamplist[[d]])),
#                                                 function(foo) matrix(gammaSamplist[[d]][foo, ], kdims[d], kdims[d]))})
#     Gkron_list <- vector("list", niter)
#     for (tt in 1:niter) {
#         obj <- lapply(seq_len(nmodes), function(d) GammaSamplist[[d]][[tt]])
#         Gkron_list[[tt]] <- revkronAll(obj)
#     }
#
#     for (gg in 1:n_test) {
#         yhat <- list()
#         for (tt in 1:niter) {
#             Gkron <- Gkron_list[[tt]]
#             di_tilde <- ranefSamplist[[tt]][gg, ] ## length should be prod(kdims)
#             ranComp <- z_tilde_list_test[[gg]] %*% Gkron %*% di_tilde
#             mu_y <- as.vector(Xi_test[[gg]] %*% betaSamp[tt, ] + ranComp)
#             yhat[[tt]] <- rnorm(mis[gg], mean = mu_y, sd = sqrt(errVarSamp[tt]))
#         }
#
#         yhatmat <- do.call(cbind, yhat)
#         if(gg == 1) {
#             y1hatmat <- yhatmat
#         }
#
#         preds_mean[[gg]] <- rowMeans(yhatmat)
#         resids_mean[[gg]] <- ylist_test[[gg]] - preds_mean[[gg]]
#
#         preds_median[[gg]] <- apply(yhatmat, 1, median)
#         resids_median[[gg]] <- ylist_test[[gg]] - preds_median[[gg]]
#
#         grpvec <- 1:mis[gg]
#
#         q025[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.025))
#         q975[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.975))
#         covgs_95[[gg]] <- sapply(grpvec, function(ii) {as.numeric(q025[[gg]][ii] <= ylist_test[[gg]][ii] & ylist_test[[gg]][ii] <= q975[[gg]][ii])})
#         ciwidth_95[[gg]] <- sapply(grpvec, function(ii) q975[[gg]][ii] - q025[[gg]][ii])
#
#         q05[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.05))
#         q95[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.95))
#         covgs_90[[gg]] <- sapply(grpvec, function(ii) {as.numeric(q05[[gg]][ii] <= ylist_test[[gg]][ii] & ylist_test[[gg]][ii] <= q95[[gg]][ii])})
#         ciwidth_90[[gg]] <- sapply(grpvec, function(ii) q95[[gg]][ii] - q05[[gg]][ii])
#
#         q10[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.10))
#         q90[[gg]] <- apply(yhatmat, 1, function(u) quantile(u, 0.90))
#         covgs_80[[gg]] <- sapply(grpvec, function(ii) {as.numeric(q10[[gg]][ii] <= ylist_test[[gg]][ii] & ylist_test[[gg]][ii] <= q90[[gg]][ii])})
#         ciwidth_80[[gg]] <- sapply(grpvec, function(ii) q90[[gg]][ii] - q10[[gg]][ii])
#
#     }
#     #preds_mean <- unlist(preds_mean); preds_median <- unlist(preds_median)
#     mspe_mean <- mean(unlist(lapply(resids_mean, function(foo) {mean(foo ^ 2)})))
#     mspe_median <- mean(unlist(lapply(resids_median, function(foo) {mean(foo ^ 2)})))
#     mad_mean <- mean(unlist(lapply(resids_mean, function(foo) {mean(abs(foo))})))
#     mad_median <- mean(unlist(lapply(resids_median, function(foo) {mean(abs(foo))})))
#
#     covgs_95 <- unlist(covgs_95)
#     ciwidth_95 <- unlist(ciwidth_95)
#
#     covgs_90 <- unlist(covgs_90)
#     ciwidth_90 <- unlist(ciwidth_90)
#
#     covgs_80 <- unlist(covgs_80)
#     ciwidth_80 <- unlist(ciwidth_80)
#
#     result <- list(preds_mean = preds_mean, mspe_mean = mspe_mean, mad_mean = mad_mean,
#                    preds_median = preds_median, mspe_median = mspe_median, mad_median = mad_median,
#                    ytrue = y_test,
#                    covgs_95 = covgs_95, ciwidth_95_summary = summary(ciwidth_95),
#                    covgs_90 = covgs_90, ciwidth_90_summary = summary(ciwidth_90),
#                    covgs_80 = covgs_80, ciwidth_80_summary = summary(ciwidth_80),
#                    q025 = unlist(q025), q975 = unlist(q975),
#                    q05 = unlist(q05), q95 = unlist(q95),
#                    q10 = unlist(q10), q90 = unlist(q90), mi_test = mis,
#                    y1hatmat = y1hatmat, R_list = R_list, S_list = S_list)
#     #y1hatmat = y1hatmat, R_list = R_list, S_list = S_list
#     result
# }

