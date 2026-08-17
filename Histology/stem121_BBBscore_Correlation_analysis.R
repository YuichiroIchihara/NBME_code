# =========================
# BBB vs STEM121 correlation
# =========================

# 必要パッケージ
library(ggplot2)

# データ入力
df <- data.frame(
  Animal = c("Purple3", "Purple4", "Purple6", "Purple8", "Black8"),
  BBB = c(7.5, 7.5, 4.0, 4.5, 2.5),
  STEM121 = c(0.441, 0.508, 0.266, 0.389, 0.256)
)

# 相関解析
pearson <- cor.test(df$BBB, df$STEM121, method = "pearson")
spearman <- cor.test(df$BBB, df$STEM121, method = "spearman", exact = FALSE)

# 結果確認
pearson
spearman

# Figure内に表示するテキスト
label_text <- paste0(
  "Pearson's r = ", round(pearson$estimate, 2),
  "\nP = ", signif(pearson$p.value, 2)
)

# Figure作成
p <- ggplot(df, aes(x = BBB, y = STEM121)) +
  geom_point(size = 3.5, color = "#2E8B57") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 2.7,
    y = 0.50,
    label = label_text,
    hjust = 0,
    vjust = 1,
    size = 4,
    family = "Arial"
  ) +
  labs(
    x = "BBB score",
    y = "STEM121-positive area"
  ) +
  theme_classic(base_family = "Arial", base_size = 14) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13, color = "black"),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    plot.title = element_blank()
  )

# 表示
p