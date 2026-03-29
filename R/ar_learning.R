run_ar3_expon_sdr <- function(spec = c("1y", "10y"),
                              year1 = 1990,
                              year2 = 2019,
                              alpha = 0.98,
                              N = 400,
                              seed = 616,
                              plot_figure = TRUE,
                              save_results = TRUE,
                              results_dir = "results",
                              fig_width = 7,
                              fig_height = 4.5) {
  
  spec <- match.arg(spec)
  set.seed(seed)
  
  plots_dir <- file.path(results_dir, "plots")
  
  if (isTRUE(save_results)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  ## load data
  if (spec == "1y") {
    data <- loadShortRate()
  } else {
    data <- loadLongRate()
  }
  
  data$r <- data$rr
  
  ## exponentially weighted moving average
  ewma <- function(x, alpha) {
    tau <- numeric(length(x))
    tau[1] <- x[1]
    for (t in 2:length(x)) {
      tau[t] <- tau[t - 1] + (1 - alpha) * (x[t] - tau[t - 1])
    }
    tau
  }
  
  mats <- 1:N
  
  ## relevant subsamples
  d1 <- data[data$year <= year1, , drop = FALSE]
  d2 <- data[data$year <= year2, , drop = FALSE]
  
  ## compute EWMA r*
  data$rstar <- ewma(data$r, alpha)
  
  rstar_year1 <- data$rstar[data$year == year1]
  rstar_year2 <- data$rstar[data$year == year2]
  
  if (length(rstar_year1) == 0) {
    stop("year1 not found in data$year: ", year1)
  }
  if (length(rstar_year2) == 0) {
    stop("year2 not found in data$year: ", year2)
  }
  
  cat("EWMA rstar", year1, ":", rstar_year1[1], "\n")
  cat("EWMA rstar", year2, ":", rstar_year2[1], "\n")
  
  ## (1) weighted means
  mm1 <- stats::lm(r ~ 1, weights = alpha^(-year), data = d1)
  etabar1 <- unname(stats::coef(mm1)[1])
  etase1 <- sqrt(sandwich::NeweyWest(mm1, prewhite = FALSE, lag = 4)[1, 1])
  cat("Mean", year1, ":", etabar1, "\n")
  cat("SE mean:", etase1, "\n")
  
  mm2 <- stats::lm(r ~ 1, weights = alpha^(-year), data = d2)
  etabar2 <- unname(stats::coef(mm2)[1])
  etase2 <- sqrt(sandwich::NeweyWest(mm2, prewhite = FALSE, lag = 4)[1, 1])
  cat("Mean", year2, ":", etabar2, "\n")
  cat("SE mean:", etase2, "\n")
  
  d1$rtilde <- d1$r - etabar1
  d2$rtilde <- d2$r - etabar2
  
  ## (2) AR(3) on residuals
  mod1 <- dynlm::dynlm(
    rtilde ~ stats::lag(rtilde, -1) + stats::lag(rtilde, -2) + stats::lag(rtilde, -3) - 1,
    data = stats::ts(d1),
    weights = tail(alpha^(year1 - d1$year), -3)
  )
  print(summary(mod1))
  cat("# Sum of AR coeffs pre-", year1, ":", sum(stats::coef(mod1)), "\n", sep = "")
  
  rhobar1 <- as.numeric(stats::coef(mod1))
  rhose1 <- sqrt(diag(sandwich::vcovHC(mod1, type = "HC0")))
  sig1 <- summary(mod1)$sigma
  
  mod2 <- dynlm::dynlm(
    rtilde ~ stats::lag(rtilde, -1) + stats::lag(rtilde, -2) + stats::lag(rtilde, -3) - 1,
    data = stats::ts(d2),
    weights = tail(alpha^(year2 - d2$year), -3)
  )
  print(summary(mod2))
  cat("# Sum of AR coeffs full sample:", sum(stats::coef(mod2)), "\n")
  
  rhobar2 <- as.numeric(stats::coef(mod2))
  rhose2 <- sqrt(diag(sandwich::vcovHC(mod2, type = "HC0")))
  sig2 <- summary(mod2)$sigma
  
  ar_tbl1 <- cbind(
    Estimate = rhobar1,
    SE = rhose1,
    `t-stat` = rhobar1 / rhose1
  )
  rownames(ar_tbl1) <- paste("lag", seq_along(rhobar1))
  print(round(ar_tbl1, 5))
  
  ar_tbl2 <- cbind(
    Estimate = rhobar2,
    SE = rhose2,
    `t-stat` = rhobar2 / rhose2
  )
  rownames(ar_tbl2) <- paste("lag", seq_along(rhobar2))
  print(round(ar_tbl2, 5))
  
  ## robust wrapper for sim_ar_sdr
  sim_formals <- names(formals(sim_ar_sdr))
  
  if (!"N" %in% sim_formals) {
    sim_env <- environment(sim_ar_sdr)
    had_N <- exists("N", envir = sim_env, inherits = FALSE)
    if (had_N) {
      old_N <- get("N", envir = sim_env, inherits = FALSE)
    }
    assign("N", N, envir = sim_env)
    
    on.exit({
      if (had_N) {
        assign("N", old_N, envir = sim_env)
      } else if (exists("N", envir = sim_env, inherits = FALSE)) {
        rm("N", envir = sim_env)
      }
    }, add = TRUE)
  }
  
  run_sim <- function(etabar, etase, rhobar, rhose, sig) {
    sim_args <- list(
      etabar = etabar,
      etase = etase,
      rhobar = rhobar,
      rhose = rhose,
      sig = sig
    )
    
    if ("N" %in% sim_formals) {
      sim_args$N <- N
    }
    
    do.call(sim_ar_sdr, sim_args)
  }
  
  ## (3) calculate term structures
  P1 <- run_sim(etabar1, etase1, rhobar1, rhose1, sig1)
  P2 <- run_sim(etabar2, etase2, rhobar2, rhose2, sig2)
  
  y1 <- -100 / mats * log(P1)
  y2 <- -100 / mats * log(P2)
  
  rstar_old <- etabar1
  rstar_new <- etabar2
  
  ## ggplot version of the same figure
  plot_df <- rbind(
    data.frame(
      Horizon = mats,
      DiscountRate = y1,
      Sample = factor(as.character(year1), levels = c(as.character(year1), as.character(year2)))
    ),
    data.frame(
      Horizon = mats,
      DiscountRate = y2,
      Sample = factor(as.character(year2), levels = c(as.character(year1), as.character(year2)))
    )
  )
  
  hline_df <- data.frame(
    Sample = factor(
      c(as.character(year1), as.character(year2)),
      levels = c(as.character(year1), as.character(year2))
    ),
    Rate = c(etabar1, etabar2)
  )
  
  cols <- stats::setNames(
    c("red", "steelblue"),
    c(as.character(year1), as.character(year2))
  )
  
  p <- ggplot(plot_df, aes(x = Horizon, y = DiscountRate, color = Sample)) +
    geom_line(linewidth = 1.1) +
    geom_hline(
      data = hline_df,
      aes(yintercept = Rate, color = Sample),
      linetype = "dashed",
      linewidth = 1.1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cols) +
    coord_cartesian(ylim = range(0, y1, y2, 4)) +
    labs(
      x = "Horizon (years)",
      y = "Discount rate (percent)",
      color = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", colour = "black"),
      legend.key = element_rect(fill = "white", colour = NA)
    )
  
  fig_file <- NULL
  data_file <- NULL
  
  if (isTRUE(plot_figure)) {
    print(p)
  }
  
  if (isTRUE(save_results)) {
    fig_file <- file.path(plots_dir, paste0("sdr_ar3_expon_", spec, ".pdf"))
    ggsave(
      filename = fig_file,
      plot = p,
      width = fig_width,
      height = fig_height,
      units = "in"
    )
    message("Figure saved to: ", fig_file)
    
    data_file <- file.path(results_dir, paste0("sdr_ar3_expon_", spec, ".RData"))
    save(
      file = data_file,
      spec, year1, year2, alpha, N, seed,
      rstar_year1, rstar_year2,
      etabar1, etase1, etabar2, etase2,
      rhobar1, rhose1, sig1, ar_tbl1,
      rhobar2, rhose2, sig2, ar_tbl2,
      rstar_old, rstar_new,
      P1, y1, P2, y2
    )
    message("Results saved to: ", data_file)
  }
  
  invisible(list(
    spec = spec,
    year1 = year1,
    year2 = year2,
    alpha = alpha,
    N = N,
    seed = seed,
    rstar_year1 = rstar_year1[1],
    rstar_year2 = rstar_year2[1],
    etabar1 = etabar1,
    etase1 = etase1,
    etabar2 = etabar2,
    etase2 = etase2,
    ar_table_pre = round(ar_tbl1, 5),
    ar_table_full = round(ar_tbl2, 5),
    rstar_old = rstar_old,
    rstar_new = rstar_new,
    P1 = P1,
    y1 = y1,
    P2 = P2,
    y2 = y2,
    plot = p,
    figure_file = fig_file,
    data_file = data_file
  ))
}