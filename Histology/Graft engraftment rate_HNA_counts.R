# ============================================================
# HNA+ cell count (engraftment) : TP only vs Combination
#  - Summary stats
#  - Mann–Whitney U test (Wilcoxon rank-sum)
#  - Bar plot (mean ± SEM) + individual dots
#    * TP left (muted orange), Combination right (green)
#    * Bars: black outline
#    * Dots: same color + black outline (shape 21)
# ============================================================

library(tidyverse)

# ----------------------------
# data
# ----------------------------
tp   <- c(13323, 12696, 8865, 9268, 7060)
comb <- c(14727, 8744, 4635, 10548, 5043)

# ----------------------------
# 1) descriptive statistics
# ----------------------------
summary_stats <- function(x){
  c(
    n      = length(x),
    mean   = mean(x),
    sd     = sd(x),
    median = median(x),
    q1     = as.numeric(quantile(x, 0.25, names = FALSE)),
    q3     = as.numeric(quantile(x, 0.75, names = FALSE))
  )
}

rbind(
  `TP only`     = summary_stats(tp),
  Combination   = summary_stats(comb)
)

# ----------------------------
# 2) Mann–Whitney U test
# ----------------------------
wilcox.test(comb, tp,
            alternative = "two.sided",
            exact = FALSE,
            correct = FALSE)

# ----------------------------
# 3) plot data frame
# ----------------------------
df <- tibble(
  group = factor(rep(c("TP only", "Combination"), each = 5),
                 levels = c("TP only", "Combination")),
  value = c(tp, comb)
)

# mean ± SEM
sum_df <- df %>%
  group_by(group) %>%
  summarise(
    mean = mean(value),
    sem  = sd(value) / sqrt(n()),
    .groups = "drop"
  )

# colors (match your style)
cols <- c(
  "TP only"     = "#D55E00",  # muted orange
  "Combination" = "#2E8B57"   # green
)

# ----------------------------
# 4) bar + dots (same design)
# ----------------------------
p <- ggplot(sum_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.15, linewidth = 0.8) +
  geom_jitter(
    data = df,
    aes(x = group, y = value, fill = group),
    width = 0.08, size = 3,
    shape = 21, color = "black", stroke = 0.8, alpha = 0.95
  ) +
  scale_fill_manual(values = cols) +
  labs(
    x = NULL,
    y = "HNA+ cell count"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

print(p)