year1 <- 1990
year2 <- 2019
N <- 400
yrange <- c(0, 4.2)

files <- c(
  "results/sdr_uc_1y.RData",
  "results/sdr_uc_10y.RData",
  "results/sdr_meanshift_1y.RData",
  "results/sdr_meanshift_10y.RData",
  "results/sdr_ar3_expon_1y.RData",
  "results/sdr_ar3_expon_10y.RData"
)

titles <- c(
  "UC model, 1y rate",
  "UC model, 10y rate",
  "AR model, break, 1y rate",
  "AR model, break, 10y rate",
  "AR model, learning, 1y rate",
  "AR model, learning, 10y rate"
)

read_panel <- function(file, title, year1 = 1990, year2 = 2019) {
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  
  line_df <- rbind(
    data.frame(
      Horizon = seq_along(e$y1),
      DiscountRate = e$y1,
      Year = factor(as.character(year1),
                    levels = c(as.character(year1), as.character(year2))),
      Panel = title
    ),
    data.frame(
      Horizon = seq_along(e$y2),
      DiscountRate = e$y2,
      Year = factor(as.character(year2),
                    levels = c(as.character(year1), as.character(year2))),
      Panel = title
    )
  )
  
  hline_df <- data.frame(
    Year = factor(c(as.character(year1), as.character(year2)),
                  levels = c(as.character(year1), as.character(year2))),
    Rate = c(e$rstar_old, e$rstar_new),
    Panel = title
  )
  
  list(lines = line_df, hlines = hline_df)
}

tmp <- Map(read_panel, files, titles)

plot_df  <- do.call(rbind, lapply(tmp, `[[`, "lines"))
hline_df <- do.call(rbind, lapply(tmp, `[[`, "hlines"))

plot_df$Panel  <- factor(plot_df$Panel, levels = titles)
hline_df$Panel <- factor(hline_df$Panel, levels = titles)

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
  facet_wrap(
    ~ Panel,
    ncol = 2,
    axes = "all",
    axis.labels = "margins"
  ) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = yrange) +
  labs(
    x = "Horizon (years)",
    y = "Discount rate (percent)",
    color = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 11),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_blank(),
    panel.spacing = grid::unit(1.1, "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.47, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_rect(fill = "white", colour = NA)
  )

plots_dir <- file.path("results", "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

fig_file <- file.path(plots_dir, "figure_2.pdf")

ggsave(
  filename = fig_file,
  plot = p,
  width = 10,
  height = 12,
  units = "in"
)

message("Figure saved to: ", fig_file)

print(p)
