sim_ar_sdr <- function(etabar, etase, rhobar, rhose, sig, M=50000, N=400) {
    df <- matrix(NA, M, N)
    for (j in 1:M) {
        rho <- rnorm(3, mean = rhobar, sd = rhose)
        while (sum(rho)>1)
            rho <- rnorm(3, mean = rhobar, sd = rhose)
        u <- numeric(N+3)
        eps <- rnorm(N, mean=0, sd=sig)
        for (t in 4:(N+3))
            u[t] = rho[1]*u[t-1] + rho[2]*u[t-2] + rho[3]*u[t-3] + eps[t-3]
        r <- rnorm(1, mean = etabar, sd = etase) + tail(u, -3)
        r[r<0] <- 0
        df[j,] <- exp(-cumsum(r/100))
    }
    colMeans(df) # discount factor
}

sim_uc_sdr <- function(rstar0, sigu, sigv, phi, LB = 0, M = 50000, N = 400) {
  ## simulate short-term SDR from UC model
  
  stopifnot(length(rstar0) == length(sigu))
  stopifnot(length(rstar0) == length(sigv))
  stopifnot(length(rstar0) == length(phi))
  stopifnot(length(rstar0) == M)
  
  rstar <- matrix(0, M, N)
  rtilde <- matrix(0, M, N)
  
  rstar[, 1] <- rstar0
  rtilde[, 1] <- 0   # no need to match current short rate
  
  for (n in 2:N) {
    rstar[, n] <- rstar[, n - 1] + rnorm(M, 0, sigu)
    rtilde[, n] <- phi * rtilde[, n - 1] + rnorm(M, 0, sigv)
  }
  
  r <- rstar + rtilde
  r[r < LB] <- LB   # shadow rate constraint
  
  discount <- t(apply(r, 1, function(x) exp(-cumsum(x / 100))))
  colMeans(discount)   # return N-vector with bond prices
}