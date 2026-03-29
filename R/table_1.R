export_sdr_table_tex <- function(year1 = 1990,
                                 year2 = 2019,
                                 results_dir = "results",
                                 out_file = file.path(results_dir, "tables", "tab_1.tex")) {
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  ## build table
  tbl <- matrix(NA_real_, 3, 6)
  rownames(tbl) <- c("UC model",
                     "AR model, break",
                     "AR model, learning")
  colnames(tbl) <- c(year1, year2, "Change",
                     year1, year2, "Change")
  
  files <- c(file.path(results_dir, "sdr_uc_1y.RData"),
             file.path(results_dir, "sdr_uc_10y.RData"),
             file.path(results_dir, "sdr_meanshift_1y.RData"),
             file.path(results_dir, "sdr_meanshift_10y.RData"),
             file.path(results_dir, "sdr_ar3_expon_1y.RData"),
             file.path(results_dir, "sdr_ar3_expon_10y.RData"))
  
  n <- 1
  for (i in 1:3) {
    for (j in 1:2) {
      e <- new.env(parent = emptyenv())
      load(files[n], envir = e)
      tbl[i, (j - 1) * 3 + (1:3)] <- c(e$rstar_old, e$rstar_new, e$rstar_new - e$rstar_old)
      n <- n + 1
    }
  }
  
  tbl <- round(tbl, 1)
  print(tbl)
  
  fmt_num <- function(x) {
    sprintf("%.1f", x)
  }
  
  tex_lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Estimates of $r_t^{*}$}",
    "\\label{tab:rstar_estimates}",
    "\\begin{threeparttable}",
    "\\setlength{\\tabcolsep}{8pt}",
    "\\begin{tabular}{lcccccc}",
    "\\toprule\\toprule",
    paste0(" & \\multicolumn{3}{c}{One-Year Rate} & \\multicolumn{3}{c}{Ten-Year Rate} \\\\"),
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    paste0("Model & ", year1, " & ", year2, " & Change & ", year1, " & ", year2, " & Change \\\\"),
    "\\midrule",
    paste0(
      rownames(tbl)[1], " & ",
      paste(fmt_num(tbl[1, ]), collapse = " & "),
      " \\\\"
    ),
    paste0(
      rownames(tbl)[2], " & ",
      paste(fmt_num(tbl[2, ]), collapse = " & "),
      " \\\\"
    ),
    paste0(
      rownames(tbl)[3], " & ",
      paste(fmt_num(tbl[3, ]), collapse = " & "),
      " \\\\"
    ),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item Model-based estimates of $r_t^{*}$ using the one-year real rate (sample: 1954--2019) or the ten-year real rate (sample: 1968--2019).",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
  
  writeLines(tex_lines, out_file)
  message("LaTeX table saved to: ", out_file)
  
  invisible(list(
    table = tbl,
    tex_file = out_file
  ))
}