# ============================================================
# Rabies-labelled reticular formation neurons
# Main Figure: raw counts
# Extended Data: normalized per 100 starter cells
# ============================================================

# グラフィックデバイスをリセット
graphics.off()

# パッケージ
library(ggplot2)

# ------------------------------------------------------------
# 1. データ入力
# ------------------------------------------------------------

df <- data.frame(
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
  RF_mCherry_cells = c(
    2, 2, 0, 0,
    0, 0, 5, 9
  ),
  Per_30_sections = c(
    2.00, 1.88, 0, 0,
    0, 0, 4.69, 8.71
  ),
  Per_100_starter_cells = c(
    0.91, 2.78, 0, 0,
    0, 0, 1.89, 3.85
  )
)

# 群の表示順を固定
df$Group <- factor(
  df$Group,
  levels = c("TP-only", "Combination")
)

# ------------------------------------------------------------
# 2. 色の指定
# ------------------------------------------------------------

group_colors <- c(
  "TP-only" = "#D55E00",
  "Combination" = "#2E8B57"
)

# ドットの重なりを避ける
jitter_position <- position_jitter(
  width = 0.08,
  height = 0,
  seed = 123
)

# 共通テーマ
common_theme <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      color = "black",
      margin = margin(r = 10)
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 13,
      margin = margin(t = 6)
    ),
    
    axis.text.y = element_text(
      size = 12
    ),
    
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
# 3. Main Figure用：網様体のrabies-labelled neuron実数
# ------------------------------------------------------------

p_main <- ggplot(
  df,
  aes(
    x = Group,
    y = RF_mCherry_cells
  )
) +
  
  # 各個体
  geom_point(
    aes(fill = Group),
    position = jitter_position,
    shape = 21,
    size = 3.5,
    color = "black",
    stroke = 0.8
  ) +
  
  # 平均 ± SEM
  stat_summary(
    aes(color = Group),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.14,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # 平均値を横線で表示
  stat_summary(
    aes(color = Group),
    fun = mean,
    geom = "point",
    shape = 95,
    size = 10,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_color_manual(
    values = group_colors
  ) +
  
  # 0のドットが切れないように、少しだけ0未満まで表示
  scale_y_continuous(
    breaks = seq(0, 10, 2),
    limits = c(-0.35, 10.5),
    expand = c(0, 0)
  ) +
  
  labs(
    y = paste0(
      "Rabies-labelled neurons\n",
      "in the reticular formation"
    )
  ) +
  
  common_theme

# Plotsペインに表示
print(p_main)

# ------------------------------------------------------------
# 4. Extended Data用：100 starter cells当たり
# ------------------------------------------------------------

p_normalized <- ggplot(
  df,
  aes(
    x = Group,
    y = Per_100_starter_cells
  )
) +
  
  # 各個体
  geom_point(
    aes(fill = Group),
    position = jitter_position,
    shape = 21,
    size = 3.5,
    color = "black",
    stroke = 0.8
  ) +
  
  # 平均 ± SEM
  stat_summary(
    aes(color = Group),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.14,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # 平均値を横線で表示
  stat_summary(
    aes(color = Group),
    fun = mean,
    geom = "point",
    shape = 95,
    size = 10,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_color_manual(
    values = group_colors
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 4, 1),
    limits = c(-0.15, 4.5),
    expand = c(0, 0)
  ) +
  
  labs(
    y = paste0(
      "Rabies-labelled neurons\n",
      "per 100 starter cells"
    )
  ) +
  
  common_theme

# Normalized graphをPlotsに表示するときに実行
# print(p_normalized)

# ------------------------------------------------------------
# 5. 統計解析
# ------------------------------------------------------------

raw_test <- wilcox.test(
  RF_mCherry_cells ~ Group,
  data = df,
  exact = FALSE
)

normalized_test <- wilcox.test(
  Per_100_starter_cells ~ Group,
  data = df,
  exact = FALSE
)

raw_test
normalized_test

# ------------------------------------------------------------
# 6. 保存
# ------------------------------------------------------------

# Main Figure
ggsave(
  filename = "RF_rabies_labelled_neurons_main.pdf",
  plot = p_main,
  width = 3.2,
  height = 4.0,
  units = "in"
)

ggsave(
  filename = "RF_rabies_labelled_neurons_main.tiff",
  plot = p_main,
  width = 3.2,
  height = 4.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# Extended Data
ggsave(
  filename = "RF_rabies_labelled_neurons_per_100_starter.pdf",
  plot = p_normalized,
  width = 3.2,
  height = 4.0,
  units = "in"
)

ggsave(
  filename = "RF_rabies_labelled_neurons_per_100_starter.tiff",
  plot = p_normalized,
  width = 3.2,
  height = 4.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

print(p_normalized)