export_scc_table_tex <- function(year1 = 1990,
                                 year2 = 2019,
                                 results_dir = "results",
                                 data_dir = "data",
                                 out_file = file.path(results_dir, "tables", "tab_2.tex")) {
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  model_files <- c(
    file.path(results_dir, "sdr_uc_1y.RData"),
    file.path(results_dir, "sdr_uc_10y.RData"),
    file.path(results_dir, "sdr_meanshift_1y.RData"),
    file.path(results_dir, "sdr_meanshift_10y.RData"),
    file.path(results_dir, "sdr_ar3_expon_1y.RData"),
    file.path(results_dir, "sdr_ar3_expon_10y.RData")
  )
  
  model_names <- c(
    "UC model, 1y rate",
    "UC model, 10y rate",
    "AR model, break, 1y rate",
    "AR model, break, 10y rate",
    "AR model, learning, 1y rate",
    "AR model, learning, 10y rate"
  )
  
  compute_panel <- function(damage_file) {
    d <- read.csv(damage_file)
    d <- d[d$year <= 2410, , drop = FALSE]
    
    mats <- d$year - d$year[1]
    loss <- d$damages
    
    stopifnot(length(loss) == length(mats))
    
    J <- length(model_files)
    tbl <- matrix(NA_real_, J, 4)
    rownames(tbl) <- model_names
    colnames(tbl) <- c("rstar_change", "SCC_1990", "SCC_2019", "pct_change")
    
    for (j in seq_len(J)) {
      e <- new.env(parent = emptyenv())
      load(model_files[j], envir = e)
      
      N <- length(e$P1)
      stopifnot(N == 400)
      stopifnot(max(mats) <= N)
      
      P1 <- c(1, e$P1)
      P2 <- c(1, e$P2)
      
      SCC1 <- as.numeric(crossprod(P1[mats + 1], loss))
      SCC2 <- as.numeric(crossprod(P2[mats + 1], loss))
      
      tbl[j, ] <- c(
        e$rstar_new - e$rstar_old,
        SCC1,
        SCC2,
        100 * (SCC2 / SCC1 - 1)
      )
    }
    
    tbl
  }
  
  tbl_h <- compute_panel(file.path(data_dir, "dice_haensel_damages.csv"))
  tbl_d <- compute_panel(file.path(data_dir, "dice_dfg_damages.csv"))
  
  print(round(tbl_h, 1))
  print(round(tbl_d, 1))
  
  fmt_num <- function(x, digits = 1) {
    formatC(round(x, digits), format = "f", digits = digits, big.mark = ",")
  }
  
  fmt_pct <- function(x) {
    paste0(formatC(round(x, 0), format = "f", digits = 0, big.mark = ","), "\\%")
  }
  
  panel_to_tex <- function(panel_title, tbl) {
    out <- c(
      paste0(panel_title, " & & & & \\\\")
    )
    
    for (i in seq_len(nrow(tbl))) {
      out <- c(
        out,
        paste0(
          "\\hspace{0.75em}", rownames(tbl)[i], " & ",
          fmt_num(tbl[i, 1], 1), " & ",
          fmt_num(tbl[i, 2], 1), " & ",
          fmt_num(tbl[i, 3], 1), " & ",
          fmt_pct(tbl[i, 4]),
          " \\\\"
        )
      )
    }
    
    out
  }
  
  tex_lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Estimates of the SCC (Dollars per Metric Ton of CO$_2$)}",
    "\\label{tab:scc_estimates}",
    "\\begin{threeparttable}",
    "\\setlength{\\tabcolsep}{6pt}",
    "\\begin{tabular}{lrrrr}",
    "\\toprule\\toprule",
    paste0("Model & Change in $r_t^{*}$ & ", year1, " & ", year2, " & \\% Change \\\\"),
    "\\midrule",
    panel_to_tex("DICE-H, from H\\\"ansel et al.\\ (2020)", tbl_h),
    "\\addlinespace[0.35em]",
    panel_to_tex("DICE-D, from Dietz et al.\\ (2020)", tbl_d),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item Estimated social cost of carbon (SCC) for six empirical SDR models and two different marginal damage profiles. The change in each model-based $r_t^{*}$ estimate from ",
    paste0(year1, " to ", year2, " is shown in percentage points. The columns ``", year1, "'' and ``", year2, "'' show the SCC using the SDR term structures for ", year1, " and ", year2, ", implied by each SDR time series model. The SCC is calculated in constant 2010 U.S. dollars, for the base year 2015, from the marginal consumption damages over 400 years resulting from one extra ton of CO$_2$ emissions. Estimated damages are based on one of two models: DICE-H (top panel) or DICE-D (lower panel)."),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
  
  writeLines(tex_lines, out_file)
  message("LaTeX table saved to: ", out_file)
  
  invisible(list(
    haensel_table = tbl_h,
    dietz_table = tbl_d,
    tex_file = out_file
  ))
}