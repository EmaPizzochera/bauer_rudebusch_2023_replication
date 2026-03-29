library(ggplot2)
run_uc_sdr <- function(spec = c("1y", "10y"),
                       seed = 616,
                       plot_figure = TRUE,
                       save_results = TRUE,
                       results_dir = "results",
                       LB = 0,
                       year1 = 1990,
                       year2 = 2019,
                       N = 400,
                       fig_width = 7,
                       fig_height = 4.5) {
  
  ## choose between 1y and 10y rates
  spec <- match.arg(spec)
  
  ## folders
  plots_dir <- file.path(results_dir, "plots")
  
  ## create output folders only if saving is requested
  if (isTRUE(save_results)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  ## load estimation results into a separate environment
  filename <- file.path(results_dir, paste0("uc_estimates_", spec, ".RData"))
  if (!file.exists(filename)) {
    stop("File not found: ", filename)
  }
  
  est <- new.env(parent = emptyenv())
  load(filename, envir = est)
  
  data <- est$data
  X <- est$X
  sigeta2_vec <- est$sigeta2_vec
  sigv2_vec <- est$sigv2_vec
  phi_vec <- est$phi_vec
  
  print(range(data$year))
  
  if (!all(c(year1, year2) %in% data$year)) {
    stop("Both year1 and year2 must be present in data$year.")
  }
  
  t1 <- which(data$year == year1)
  t2 <- which(data$year == year2)
  
  tbl <- cbind(data$rr[c(t1, t2)], data$rstar.mean[c(t1, t2)])
  tbl <- cbind(tbl, tbl[, 1] - tbl[, 2])
  colnames(tbl) <- c("Real rate", "r*", "r-tilde")
  rownames(tbl) <- c(year1, year2)
  print(round(tbl, 2))
  
  rstar_old <- mean(X[, 1, t1])
  rstar_new <- mean(X[, 1, t2])
  
  ## term structure of discount rates
  mats <- 1:N
  cat("Simulating from posterior distribution...\n")
  
  set.seed(seed)
  P1 <- sim_uc_sdr(
    X[, 1, t1],
    sqrt(sigeta2_vec),
    sqrt(sigv2_vec),
    phi_vec,
    LB = LB,
    M = length(X[, 1, t1]),
    N = N
  )
  
  P2 <- sim_uc_sdr(
    X[, 1, t2],
    sqrt(sigeta2_vec),
    sqrt(sigv2_vec),
    phi_vec,
    LB = LB,
    M = length(X[, 1, t2]),
    N = N
  )
  
  y1 <- -100 / (1:N) * log(P1)
  f1 <- -100 * diff(log(c(1, P1)))
  y2 <- -100 / (1:N) * log(P2)
  f2 <- -100 * diff(log(c(1, P2)))
  
  ## data for ggplot
  plot_df <- rbind(
    data.frame(
      Horizon = mats,
      DiscountRate = y1,
      Year = factor(as.character(year1), levels = c(as.character(year1), as.character(year2)))
    ),
    data.frame(
      Horizon = mats,
      DiscountRate = y2,
      Year = factor(as.character(year2), levels = c(as.character(year1), as.character(year2)))
    )
  )
  
  hline_df <- data.frame(
    Year = factor(
      c(as.character(year1), as.character(year2)),
      levels = c(as.character(year1), as.character(year2))
    ),
    Rate = c(rstar_old, rstar_new)
  )
  
  cols <- stats::setNames(
    c("red", "steelblue"),
    c(as.character(year1), as.character(year2))
  )
  
  p <- ggplot(plot_df, aes(x = Horizon, y = DiscountRate, color = Year)) +
    geom_line(linewidth = 1.1) +
    geom_hline(
      data = hline_df,
      aes(yintercept = Rate, color = Year),
      linetype = "dashed",
      linewidth = 1.1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cols) +
    coord_cartesian(ylim = range(0, y1, y2)) +
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
  
  ## display figure
  if (isTRUE(plot_figure)) {
    print(p)
  }
  
  ## save figure with ggsave
  if (isTRUE(save_results)) {
    fig_file <- file.path(plots_dir, paste0("sdr_uc_", spec, ".pdf"))
    ggsave(
      filename = fig_file,
      plot = p,
      width = fig_width,
      height = fig_height,
      units = "in"
    )
    message("Figure saved to: ", fig_file)
  }
  
  ## save term structures only if requested
  if (isTRUE(save_results)) {
    data_file <- file.path(results_dir, paste0("sdr_uc_", spec, ".RData"))
    save(
      file = data_file,
      rstar_old, rstar_new, P1, f1, y1, P2, f2, y2, N, seed
    )
    message("Results saved to: ", data_file)
  }
  
  invisible(list(
    spec = spec,
    seed = seed,
    table = round(tbl, 2),
    rstar_old = rstar_old,
    rstar_new = rstar_new,
    P1 = P1, f1 = f1, y1 = y1,
    P2 = P2, f2 = f2, y2 = y2,
    N = N,
    plot = p,
    figure_file = fig_file,
    data_file = data_file
  ))
}