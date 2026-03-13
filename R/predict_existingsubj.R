#' @title predict_existingsubj
#'
#' @description This function performs posterior prediction of new observations for existing subjects using the CoMET model.
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

predict_existingsubj <- function(object, kdims,
                                 R_list, S_list,
                                 xlist_test, zlist_test, mis,
                                 nom.level) {

    betaSamp <- object$betaSamp
    errVarSamp <- object$errVarSamp
    gammaSamplist <- object$gammaSamplist
    ranefSamplist <- object$ranefSamplist

    nmodes <- length(gammaSamplist)
    n_test <- length(mis); N_test <- sum(mis)
    mis_cumsum <- cumsum(mis)
    mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
    niter <- length(errVarSamp)

    resids <- vector("list", n_test); preds <- vector("list", n_test)
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

    vecx_test <- lapply(xlist_test, function(foo) {as.vector(foo)})
    Xi_test <- lapply(1:n_test,
                      function(i) {
                          rows_i <- mis_starts[i]:mis_cumsum[i]
                          do.call("rbind", vecx_test[rows_i])
                      })

    GammaSamplist <- lapply(seq_len(nmodes),
                            function(d) {lapply(seq_len(nrow(gammaSamplist[[d]])),
                                                function(foo) matrix(gammaSamplist[[d]][foo, ], kdims[d], kdims[d]))})
    Gkron_list <- vector("list", niter)
    for (tt in 1:niter) {
        obj <- lapply(seq_len(nmodes), function(d) GammaSamplist[[d]][[tt]])
        Gkron_list[[tt]] <- revkronAll(obj)
    }

    for (gg in 1:n_test) {
        yhat <- list()
        for (tt in 1:niter) {
            Gkron <- Gkron_list[[tt]]
            di_tilde <- ranefSamplist[[tt]][[gg]]
            ranComp <- z_tilde_list_test[[gg]] %*% Gkron %*% di_tilde
            mu_y <- as.vector(Xi_test[[gg]] %*% betaSamp[tt, ] + ranComp)
            yhat[[tt]] <- rnorm(mis[gg], mean = mu_y, sd = sqrt(errVarSamp[tt]))
        }

        yhat_samples <- do.call(cbind, yhat)

        preds[[gg]] <- rowMeans(yhat_samples)

        grpvec <- 1:mis[gg]

        qlower[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, (1 - nom.level) / 2))
        qupper[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, 1 - (1 - nom.level) / 2))

    }

    result <- list(ypred = preds, yhat_samples = yhat_samples,
                   lower_pi = qlower, upper_pi = qupper)

    result
}

