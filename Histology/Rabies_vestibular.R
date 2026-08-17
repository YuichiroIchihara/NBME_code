# ============================================================
# Vestibular nuclei rabies-labelled neurons
# Raw count and normalized per 100 starter cells
# ============================================================

graphics.off()
library(ggplot2)

# ------------------------------------------------------------
# 1. データ入力
# ------------------------------------------------------------

df_vestibular <- data.frame(
  Group = c(
    "TP-only", "TP-only", "TP-only", "TP-only",
    "Combination", "Combination", "Combination", "Combination"
  ),
  Animal = c(
    "G1", "G5", "G7", "G6",
    "G11", "G4", "G12", "G8"
  ),
  Starter_cells = c(
    219, 72, 58, 53,
    98, 51, 265, 234
  ),
  Analyzed_sections = c(
    30, 32, 29, 32,
    32, 32, 32, 31
  ),
  Vestibular_mCherry_cells = c(
    0, 1, 0, 0,
    3, 0, 0, 0
  )
)

# 100 starter cells当たり
df_vestibular$Vestibular_per_100_starter <-
  df_vestibular$Vestibular_mCherry_cells /
  df_vestibular$Starter_cells * 100

# 群の順番
df_vestibular$Group <- factor(
  df_vestibular$Group,
  levels = c("TP-only", "Combination")
)

# ------------------------------------------------------------
# 2. 色
# ------------------------------------------------------------

group_colors <- c(
  "TP-only" = "#D55E00",
  "Combination" = "#2E8B57"
)

jitter_position <- position_jitter(
  width = 0.08,
  height = 0,
  seed = 123
)

# ------------------------------------------------------------
# 3. 共通テーマ
# ------------------------------------------------------------

common_theme <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      color = "black",
      margin = margin(r = 10)
    ),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(
      size = 13,
      margin = margin(t = 6)
    ),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.7
    ),
    axis.ticks.length = unit(2.5, "mm"),
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 10,
      l = 12
    )
  )

# ------------------------------------------------------------
# 4. 前庭核 raw count
# ------------------------------------------------------------

p_vestibular_raw <- ggplot(
  df_vestibular,
  aes(
    x = Group,
    y = Vestibular_mCherry_cells
  )
) +
  geom_point(
    aes(fill = Group),
    position = jitter_position,
    shape = 21,
    size = 3.5,
    color = "black",
    stroke = 0.8
  ) +
  stat_summary(
    aes(color = Group),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.14,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  stat_summary(
    aes(color = Group),
    fun = mean,
    geom = "point",
    shape = 95,
    size = 10,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    breaks = seq(0, 4, 1),
    limits = c(-0.15, 4),
    expand = c(0, 0)
  ) +
  labs(
    y = paste0(
      "Rabies-labelled neurons\n",
      "in the vestibular nuclei"
    )
  ) +
  common_theme

# Plotsに表示
print(p_vestibular_raw)

# ------------------------------------------------------------
# 5. 前庭核 per 100 starter cells
# ------------------------------------------------------------

p_vestibular_norm <- ggplot(
  df_vestibular,
  aes(
    x = Group,
    y = Vestibular_per_100_starter
  )
) +
  geom_point(
    aes(fill = Group),
    position = jitter_position,
    shape = 21,
    size = 3.5,
    color = "black",
    stroke = 0.8
  ) +
  stat_summary(
    aes(color = Group),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.14,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  stat_summary(
    aes(color = Group),
    fun = mean,
    geom = "point",
    shape = 95,
    size = 10,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    breaks = seq(0, 4, 1),
    limits = c(-0.15, 4),
    expand = c(0, 0)
  ) +
  labs(
    y = paste0(
      "Rabies-labelled neurons\n",
      "per 100 starter cells"
    )
  ) +
  common_theme

# 正規化グラフをPlotsに表示するとき
# print(p_vestibular_norm)

# ------------------------------------------------------------
# 6. 統計
# ------------------------------------------------------------

vestibular_raw_test <- wilcox.test(
  Vestibular_mCherry_cells ~ Group,
  data = df_vestibular,
  exact = FALSE
)

vestibular_norm_test <- wilcox.test(
  Vestibular_per_100_starter ~ Group,
  data = df_vestibular,
  exact = FALSE
)

vestibular_raw_test
vestibular_norm_test

# ------------------------------------------------------------
# 7. 保存
# ------------------------------------------------------------

ggsave(
  filename = "Vestibular_nuclei_rabies_labelled_neurons_raw.pdf",
  plot = p_vestibular_raw,
  width = 3.2,
  height = 4.0,
  units = "in"
)

ggsave(
  filename = "Vestibular_nuclei_rabies_labelled_neurons_raw.tiff",
  plot = p_vestibular_raw,
  width = 3.2,
  height = 4.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Vestibular_nuclei_rabies_labelled_neurons_per_100_starter.pdf",
  plot = p_vestibular_norm,
  width = 3.2,
  height = 4.0,
  units = "in"
)

ggsave(
  filename = "Vestibular_nuclei_rabies_labelled_neurons_per_100_starter.tiff",
  plot = p_vestibular_norm,
  width = 3.2,
  height = 4.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)