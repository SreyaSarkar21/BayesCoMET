#' @title predict_newsubj
#'
#' @description This function performs posterior prediction for new subjects using the compressed mixed model.
#'
#' @param object a list of posterior samples of parameters returned by the function \code{\link{comet}}.
#' @param kdims a vector of length \eqn{D} (number of tensor modes), where each element equals the \eqn{d}-th mode-specific compressed covariance dimension \eqn{k_d}.
#' @param R_list the list of \eqn{k_d \times q} random projection matrices used for compressing the random slope in the Gibbs sampler of CoMET. The list should be same as that used as argument to the function \code{\link{comet}}.
#' @param S_list the list of \eqn{k_d \times q} random projection matrices used for compressing the random-effect covariates in the Gibbs sampler of CoMET. The list should be same as that used as argument to the function \code{\link{comet}}.
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
                            xlist_test, zlist_test, mis, nom.level) {

    nmodes <- length(kdims)
    if(nmodes == 1){
        xlist_test_cme <- lapply(split(xlist_test, rep(seq_along(mis), times = mis)), function(foo) do.call("rbind", foo))
        xlist_test_cme <- unname(xlist_test_cme)
        zlist_test_cme <- lapply(split(zlist_test, rep(seq_along(mis), times = mis)), function(foo) do.call("rbind", foo))
        zlist_test_cme <- unname(zlist_test_cme)

        result <- predictCME_newsubj(sampler_res = object,
                                       xlist_test = xlist_test_cme, zlist_test = zlist_test_cme,
                                       R = R_list[[1]], S = S_list[[1]], nom.level = nom.level)

    } else{

    betaSamp <- object$betaSamp
    errVarSamp <- object$errVarSamp
    gammaSamp <- object$gammaSamp

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
                            function(d) {lapply(seq_len(nrow(gammaSamp[[d]])),
                                                function(foo) matrix(gammaSamp[[d]][foo, ], kdims[d], kdims[d]))})

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
    }

    result
}


