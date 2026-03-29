run_ar_meanshift_sdr <- function(spec = c("1y", "10y"),
                                 byear = 1990,
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
  
  ## load data depending on selected spec
  if (spec == "1y") {
    data <- loadShortRate()
  } else {
    data <- loadLongRate()
  }
  
  mats <- 1:N
  
  ## subtract out different unconditional mean pre and post break
  d1 <- data[data$year <= byear, , drop = FALSE]
  d2 <- data[data$year > byear, , drop = FALSE]
  
  mean1 <- mean(d1$rr, na.rm = TRUE)
  mean2 <- mean(d2$rr, na.rm = TRUE)
  
  data$rtilde <- data$rr - ifelse(data$year <= byear, mean1, mean2)
  
  mm <- stats::lm(rr ~ 1, data = d1)
  etabar1 <- unname(stats::coef(mm)[1])
  etase1 <- sqrt(sandwich::NeweyWest(mm)[1, 1])
  cat("Mean and SE before break:", etabar1, etase1, "\n")
  
  mm <- stats::lm(rr ~ 1, data = d2)
  etabar2 <- unname(stats::coef(mm)[1])
  etase2 <- sqrt(sandwich::NeweyWest(mm)[1, 1])
  cat("Mean and SE after break:", etabar2, etase2, "\n")
  
  stopifnot(all.equal(mean1, etabar1, check.attributes = FALSE))
  stopifnot(all.equal(mean2, etabar2, check.attributes = FALSE))
  
  mm <- stats::lm(rr ~ I(year <= byear), data = data)
  stopifnot(all.equal(mean2, unname(stats::coef(mm)[1]), check.attributes = FALSE))
  
  mean_change_tstat <- unname(stats::coef(mm)[2] / sqrt(diag(sandwich::NeweyWest(mm)))[2])
  cat("t-stat on mean change:", mean_change_tstat, "\n")
  
  ## AR model on residuals
  mod <- stats::ar.ols(data$rtilde, aic = FALSE, order.max = 3, demean = FALSE)
  rhobar <- as.numeric(drop(mod$ar))
  rhose <- as.numeric(drop(mod$asy.se.coef$ar))
  sig <- sqrt(mod$var.pred)
  
  cat("# Sum of AR coeffs:", sum(rhobar), "\n")
  
  ar_tbl <- cbind(
    Estimate = rhobar,
    SE = rhose,
    `t-stat` = rhobar / rhose
  )
  rownames(ar_tbl) <- paste("lag", seq_along(rhobar))
  print(round(ar_tbl, 5))
  
  ## simulate SDR
  sim_formals <- names(formals(sim_ar_sdr))
  
  ## this makes the wrapper work whether sim_ar_sdr() already accepts N
  ## or still relies on a global N
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
  
  run_sim <- function(mu, mu_se) {
    sim_args <- list(mu, mu_se, rhobar, rhose, sig)
    if ("N" %in% sim_formals) {
      sim_args$N <- N
    }
    do.call(sim_ar_sdr, sim_args)
  }
  
  P1 <- run_sim(etabar1, etase1)
  P2 <- run_sim(etabar2, etase2)
  
  y1 <- -100 / mats * log(P1)
  y2 <- -100 / mats * log(P2)
  
  rstar_old <- mean1
  rstar_new <- mean2
  
  ## ggplot version of the same figure
  pre_label <- paste0("pre-", byear, " mean")
  post_label <- paste0("post-", byear, " mean")
  
  plot_df <- rbind(
    data.frame(
      Horizon = mats,
      DiscountRate = y1,
      Regime = factor(pre_label, levels = c(pre_label, post_label))
    ),
    data.frame(
      Horizon = mats,
      DiscountRate = y2,
      Regime = factor(post_label, levels = c(pre_label, post_label))
    )
  )
  
  hline_df <- data.frame(
    Regime = factor(c(pre_label, post_label), levels = c(pre_label, post_label)),
    Rate = c(mean1, mean2)
  )
  
  cols <- stats::setNames(
    c("red", "steelblue"),
    c(pre_label, post_label)
  )
  
  p <- ggplot(plot_df, aes(x = Horizon, y = DiscountRate, color = Regime)) +
    geom_line(linewidth = 1.1) +
    geom_hline(
      data = hline_df,
      aes(yintercept = Rate, color = Regime),
      linetype = "dashed",
      linewidth = 1.1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cols) +
    coord_cartesian(ylim = range(0, y1, y2, 5)) +
    labs(
      x = "Horizon (years)",
      y = "Discount rate (percent",
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
    fig_file <- file.path(plots_dir, paste0("sdr_meanshift_", spec, ".pdf"))
    ggsave(
      filename = fig_file,
      plot = p,
      width = fig_width,
      height = fig_height,
      units = "in"
    )
    message("Figure saved to: ", fig_file)
    
    data_file <- file.path(results_dir, paste0("sdr_meanshift_", spec, ".RData"))
    save(
      file = data_file,
      byear, N, seed,
      mean1, mean2,
      etabar1, etase1, etabar2, etase2,
      mean_change_tstat,
      rhobar, rhose, sig, ar_tbl,
      rstar_old, rstar_new,
      P1, y1, P2, y2
    )
    message("Results saved to: ", data_file)
  }
  
  invisible(list(
    spec = spec,
    byear = byear,
    N = N,
    seed = seed,
    mean1 = mean1,
    mean2 = mean2,
    mean_change_tstat = mean_change_tstat,
    ar_table = round(ar_tbl, 5),
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