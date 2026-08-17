library(tidyverse)

# data
df <- tibble(
  group = rep(c("Control", "EESREHAB only"), each = 4),
  ratio = c(
    0.75862069, 0.742424242, 0.690217391, 0.682795699,
    1.017857143, 1.054945055, 1.064705882, 1.069364162
  )
) %>%
  mutate(group = factor(group, levels = c("Control", "EESREHAB only")))

# Welch t-test (two-tailed)
tt <- t.test(ratio ~ group, data = df, var.equal = FALSE)
pval <- tt$p.value
print(tt)

# summary: mean ± SEM
sum_df <- df %>%
  group_by(group) %>%
  summarise(
    mean = mean(ratio),
    sd   = sd(ratio),
    n    = n(),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

# colors
cols <- c("Control" = "#4D4D4D", "EESREHAB only" = "#E377C2")

# y headroom
ymax <- max(df$ratio) * 1.12

# plot
p <- ggplot(sum_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.60, color = "black", linewidth = 0.7) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.18, linewidth = 0.7
  ) +
  geom_point(
    data = df,
    aes(x = group, y = ratio, fill = group),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.10, height = 0),
    shape = 21, color = "black", stroke = 0.7,
    size = 2.6, alpha = 0.95
  ) +
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(0, ymax)) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.line    = element_line(linewidth = 0.9),
    axis.ticks   = element_line(linewidth = 0.9)
  ) +
  # significance bracket + "*"
  annotate("segment", x = 1, xend = 2, y = ymax*0.97, yend = ymax*0.97, linewidth = 0.9) +
  annotate("segment", x = 1, xend = 1, y = ymax*0.95, yend = ymax*0.97, linewidth = 0.9) +
  annotate("segment", x = 2, xend = 2, y = ymax*0.95, yend = ymax*0.97, linewidth = 0.9) +
  annotate("text", x = 1.5, y = ymax*1.01,
           label = ifelse(pval < 0.05, "*", "n.s."),
           family = "Arial", size = 6)

p
