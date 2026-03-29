get_prior <- function(p=1) {
    prior <- list(x0 = c(2, rep(0, p)), V0 = 100*diag(1+p))
    ## proper priors for variances:
    ## IG(shape, rate), shape = alpha/2, rate = delta/2, mean = delta/(alpha-2), mode = delta/(alpha+2)
    prior <- within(prior, {
        ## phi ~ TN(0, var_phi) truncated to [-1,1]
        var_phi = 1
        ## r* innovation variance -- tight prior
        ## mean matches mode of Del Negro
        ## (they use prior mode 0.0025 for r* and 0.01 for pi* innovation var. in quarterly data)
        alpha_eta = 100
        delta_eta = 0.04 * (alpha_eta - 2)  # note: annual as opposed to quarterly
        ## AR(1) innovation variance
        mean_sigv2 = 1
        alpha_v = 6
        delta_v = mean_sigv2*(alpha_v - 2)})
    prior$mean_sigv2 <- NULL
    prior
}

getRandomPars <- function(p = 1) {
    list(phi = runif(1, -1, 1),
         sigv2 = 1/rgamma(1, 1),
         sigeta2 = 1/rgamma(1, 1))
}

stateSpaceMats <- function(pars, a0, P0) {
    ## uCovVar <- function(phi, omega) {
    ##     ## unconditional covariance matrix of VAR(1)
    ##     ## - can also be used to find unconditional variance of AR(p) process by using companion form
    ##     N <- nrow(omega)
    ##     matrix(solve(diag(N^2) - kronecker(phi, phi)) %*% as.numeric(omega), N, N)
    ## }
    p <- length(pars$phi)
    N <- p + 1
    if (length(pars$theta)>0)
        stop("this is for AR(p) only, no MA terms allowed")
    ## measurement equation
    ct <- matrix(0, 1, 1)
    Z <- matrix(c(1, 1, rep(0, p-1)), 1, N)
    Zt <- array(Z, c(1, N, 1))
    GGt <- array(0, c(1, 1, 1))
    ## transition equation
    dt <- matrix(0, N, 1)
    if (p==1) {
        T <- diag(c(1, pars$phi))
    } else {
        T <- rbind(c(1, rep(0, p)),
                   c(0, pars$phi),
                   cbind(0, diag(p-1), 0))
    }
    Tt <- array(T, c(N, N, 1))
    R <- rbind(diag(2), matrix(0, p-1, 2))
    Q <- diag(c(pars$sigeta2, pars$sigv2))
    HHt <- array(R %*% Q %*% t(R), c(N, N, 1))
    if (missing(a0))
        a0 <- c(2, rep(0, p))
    if (missing(P0))
        P0 <- 100*diag(1+p)
    list(a0=a0, P0=P0, ct=ct, Z=Z, Zt=Zt, H=0, GGt=GGt, dt=dt, T=T, Tt=Tt, R=R, Q=Q, HHt=HHt)
}

ssm_uc <- function(y, pars, prior) {
    ## state-space model - UC model
    ss <- stateSpaceMats(pars, prior$x0, prior$V0)
    ## state-space model
    ## Measurement equation:
    ## y_t = Z * alpha_t + w_t              E(ww') = H
    ## Transition equation:
    ## alpha_t = T * alpha_t-1 + R*u_t      E(uu') = Q
    SSModel(y ~ -1 + SSMcustom(Z = ss$Z, T = ss$T, R = ss$R, Q = ss$Q, a1 = ss$a0, P1 = ss$P0), H = ss$H)
}

draw_regression <- function(ydat, Xdat, sig2, prior_mean, prior_precision, restrict) {
    ## draw from conditional posterior (given the data and sig2)
    ## - prior for coefficients is N(beta0, sigma^2*B0)
    ## - prior_precision = B0^-1
    ## - flag restrict=TRUE -> make sure last coefficient is between -1 and 1, so that AR(1) is stationary
    B1 <- solve(1/sig2 * t(Xdat) %*% Xdat + prior_precision)
    b_post <- B1 %*% (1/sig2 * t(Xdat) %*% ydat + prior_precision %*% prior_mean)
    valid <- FALSE
    while (!valid) {
        res <- drawNormal(b_post, B1)
        valid <- ifelse(restrict, abs(tail(res,1))<1, TRUE)
    }
    drop(res)
}

drawNormal <- function(mu, Omega)
    mu + t(chol(Omega)) %*% rnorm(length(mu))
