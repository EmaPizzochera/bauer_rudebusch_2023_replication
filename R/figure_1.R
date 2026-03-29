## plot r-star estimates from UC model
library(dplyr)
library(ggplot2)
library(patchwork)
q <- 0.68   ## 68 percent credibility intervals (1 SE width for normal distribution)

cat(100*q, "percent CIs\n")
cat("Lower quantile:", lq <- (1-q)/2, "\n")
cat("Upper quantile:", uq <- (1+q)/2, "\n")

load("results/uc_estimates_1y.RData")
stopifnot(all.equal(data$rstar.mean, colMeans(X[,1,])))
data <- data %>%
    rename(r1 = rr,
           rstar.r1 = rstar.mean) %>%
    mutate(rstar.r1.lb = apply(X[,1,], 2, quantile, lq),
           rstar.r1.ub = apply(X[,1,], 2, quantile, uq)) %>%
    select(year, r1, starts_with("rstar.r1"))
env <- new.env()
load("results/uc_estimates_10y.RData", envir=env)
d10 <- env$data %>%
    rename(r10 = rr,
           rstar.r10 = rstar.mean) %>%
    mutate(rstar.r10.lb = apply(env$X[,1,], 2, quantile, lq),
           rstar.r10.ub = apply(env$X[,1,], 2, quantile, uq)) %>%
    select(year, r10, starts_with("rstar.r10"))
data <- left_join(data, d10)
data %>% filter(!is.na(r1)) %>% pull(year) %>% range
data %>% filter(!is.na(r10)) %>% pull(year) %>% range

## external r-star estimates
data <- data %>%
    mutate(yyyymm = year*100+12) %>%
    left_join(loadDN()) %>%
    left_join(loadJM()) %>%
    left_join(loadLW()) %>%
    left_join(loadKiley()) %>%
    select(-yyyymm)
nms <- c("rstar.dn.sm", "rstar.jm.sm", "rstar.lw.sm", "rstar.kiley.sm")
data$rstar_mean <- rowMeans(as.matrix(data[nms]))

## Figure 1 - rstar from UC models for 1y and 10y rate
dev.new()
## common y-range, same idea as base R
yrange <- data %>%
  select(-year) %>%
  range(na.rm = TRUE)

## left panel: 1-year real rate
p1 <- ggplot(data, aes(x = year)) +
  geom_line(aes(y = r1, color = "Real rate"), linewidth = 0.5) +
  geom_line(aes(y = rstar.r1, color = "UC model r*"), linewidth = 1.4) +
  geom_line(aes(y = rstar.r1.lb), linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = rstar.r1.ub), linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = rstar_mean, color = "Avg. other r*"), linewidth = 0.9, linetype = "longdash") +
  scale_color_manual(
    values = c("Real rate" = "gray", "UC model r*" = "black", "Avg. other r*" = "red")
  ) +
  coord_cartesian(xlim = range(data$year, na.rm = TRUE), ylim = yrange, expand = FALSE) +
  labs(
    x = "Year",
    y = "Percent",
    title = "One-year real rate",
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "plain", size = 11),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

## right panel: 10-year real rate
d2 <- filter(data, !is.na(r10))

p2 <- ggplot(d2, aes(x = year)) +
  geom_line(aes(y = r10), color = "gray", linewidth = 0.5) +
  geom_line(aes(y = rstar.r10), color = "black", linewidth = 1.4) +
  geom_line(aes(y = rstar.r10.lb), color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = rstar.r10.ub), color = "black", linetype = "dashed", linewidth = 0.5) +
  coord_cartesian(xlim = range(d2$year, na.rm = TRUE), ylim = yrange, expand = FALSE) +
  labs(
    x = "Year",
    y = NULL,
    title = "Ten-year real rate"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "plain", size = 11),
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 6))
  )

## combine side by side
p <- p1 + p2 + plot_layout(ncol = 2)

#Save automatically in pdf format
ggsave("results/plots/figure_1.pdf", plot = p, width = 10, height = 4)
