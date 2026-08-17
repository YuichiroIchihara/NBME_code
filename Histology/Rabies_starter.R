# ============================================================
# Starter cell count
# Individual animals with mean ± SEM
# ============================================================

graphics.off()
library(ggplot2)

# データ
df_starter <- data.frame(
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
  )
)

# 群の順番
df_starter$Group <- factor(
  df_starter$Group,
  levels = c("TP-only", "Combination")
)

# 群の色
group_colors <- c(
  "TP-only" = "#D55E00",
  "Combination" = "#2E8B57"
)

# グラフ
p_starter <- ggplot(
  df_starter,
  aes(x = Group, y = Starter_cells)
) +
  
  # 各個体
  geom_point(
    aes(fill = Group),
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 123
    ),
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
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  scale_y_continuous(
    breaks = seq(0, 300, 50),
    limits = c(0, 300),
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  labs(
    x = NULL,
    y = "Starter cells"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    
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

# Plotsペインに表示
print(p_starter)

starter_test <- wilcox.test(
  Starter_cells ~ Group,
  data = df_starter,
  exact = FALSE
)

starter_test