# Mann–Whitney U test (Wilcoxon rank-sum) for 2 groups
# Vsx2周囲ROI内シナプス数：Combination vs TP

comb <- c(10.34, 17.34, 8.00, 1.34, 11.34)
tp   <- c(4.33, 0.00, 2.67, 0.34, 0.00)

# 1) 記述統計（中央値・IQR など）
summary_stats <- function(x){
  c(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    q1 = quantile(x, 0.25, names = FALSE),
    q3 = quantile(x, 0.75, names = FALSE)
  )
}
rbind(Combination = summary_stats(comb),
      TP          = summary_stats(tp))

# 2) Mann–Whitney U（Wilcoxon rank-sum）
# exact=FALSE：同順位(ties)があり得る/小数でexact計算が不安定な場合にも安定
# correct=FALSE：連続補正なし（入れたいならTRUE）
wilcox.test(comb, tp,
            alternative = "two.sided",
            exact = FALSE,
            correct = FALSE)


# ----------------------------
# Bar plot + dots (TP left / Combination right)
# TP: muted orange, Combination: green
# Bars: black outline
# Dots: same color + black outline
# ----------------------------

library(tidyverse)

tp   <- c(4.33, 0.00, 2.67, 0.34, 0.00)
comb <- c(10.34, 17.34, 8.00, 1.34, 11.34)

df <- tibble(
  group = factor(rep(c("TP", "Combination"), each = 5),
                 levels = c("TP", "Combination")),
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

# 色（渋めオレンジ & 緑）
cols <- c(
  "TP" = "#D55E00",          # muted orange（渋め）
  "Combination" = "#2E8B57"  # green
)

ggplot(sum_df, aes(x = group, y = mean, fill = group)) +
  # 棒グラフ：黒縁つき
  geom_col(width = 0.6, color = "black", linewidth = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.15, linewidth = 0.8) +
  # ドット：同色 + 黒縁
  geom_jitter(data = df, aes(x = group, y = value, fill = group),
              width = 0.08, size = 3,
              shape = 21, color = "black", stroke = 0.8, alpha = 0.95) +
  scale_fill_manual(values = cols) +
  labs(
    x = NULL,
    y = "Synaptic puncta count (Vsx2 ROI)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")