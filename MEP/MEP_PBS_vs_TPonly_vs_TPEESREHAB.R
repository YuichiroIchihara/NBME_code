# analysis title
# date

#Clear R's brain
rm(list=ls()) 

# データ入力
pbs  <- c(6, 11, 13)
comb <- c(65, 52, 51)

df <- data.frame(
  group = factor(rep(c("PBS", "Combination"), each = 3)),
  value = c(pbs, comb)
)

# 必要パッケージ
library(ggplot2)
library(dplyr)

# データ入力
pbs  <- c(6, 11, 13, 40, 31, 53)
tp   <- c(11, 19, 56, 43, 53)
comb <- c(65, 52, 51, 146, 84, 116)

df <- data.frame(
  group = factor(c(rep("PBS", 6),
                   rep("TP", 5),
                   rep("Combination", 6)),
                 levels = c("PBS", "TP", "Combination")),
  value = c(pbs, tp, comb)
)

# 色の設定（PBS：ピンク / TP：グレー / Combination：シアン）
pbs_color  <- "#0072B2"
tp_color   <- "#D55E00"
comb_color <- "#2E8B57"

# mean + SEM
summary_df <- df %>%
  group_by(group) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value),
    n    = length(value),
    sem  = sd / sqrt(n)
  )

# プロット
ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = mean - sem,
                    ymax = mean + sem),
                width = 0.2, linewidth = 1) +
  geom_jitter(data = df,
              aes(x = group, y = value, fill=group),
              width = 0.08,
              size = 3,
              shape = 21,
              color = "black",
              stroke = 1.2) +
  scale_fill_manual(values = c("PBS" = pbs_color,
                               "TP" = tp_color,
                               "Combination" = comb_color)) +
  labs(y = "MEP amplitude (mV)", x = "") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14)
  )

## =========================
## 統計解析
## =========================

# 必要パッケージ
if (!requireNamespace("rstatix", quietly = TRUE)) {
  install.packages("rstatix")
}
library(rstatix)

## -------------------------
## 1) 記述統計
## -------------------------
summary_df
# または
df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    sem = sd / sqrt(n),
    median = median(value),
    IQR = IQR(value),
    .groups = "drop"
  )

## -------------------------
## 2) 正規性の確認
## 各群 n が小さいので参考程度
## -------------------------
shapiro_res <- df %>%
  group_by(group) %>%
  shapiro_test(value)

shapiro_res

## -------------------------
## 3) 等分散性の確認
## -------------------------
levene_res <- df %>%
  levene_test(value ~ group)

levene_res

## =========================
## パラメトリック解析
## =========================

## -------------------------
## 4A) 通常の one-way ANOVA
## 分散が同程度とみなす場合
## -------------------------
anova_res <- aov(value ~ group, data = df)
summary(anova_res)

## Tukey post hoc
tukey_res <- TukeyHSD(anova_res)
tukey_res


## -------------------------
## 4B) Welch ANOVA
## 分散が異なる可能性を考慮
## 個人的にはこのデータではこちらを優先
## -------------------------
welch_res <- df %>%
  welch_anova_test(value ~ group)

welch_res

## Games-Howell post hoc
## Welch ANOVA後の多重比較として自然
games_res <- df %>%
  games_howell_test(value ~ group)

games_res


## -------------------------
## 4C) 参考：Welch t-testのpairwise比較
## Holm補正
## -------------------------
pairwise_welch_res <- pairwise.t.test(
  x = df$value,
  g = df$group,
  p.adjust.method = "holm",
  pool.sd = FALSE
)

pairwise_welch_res


## =========================
## ノンパラメトリック解析
## =========================

## -------------------------
## 5A) Kruskal-Wallis test
## -------------------------
kruskal_res <- df %>%
  kruskal_test(value ~ group)

kruskal_res

## -------------------------
## 5B) Dunn test
## Holm補正
## -------------------------
dunn_res <- df %>%
  dunn_test(value ~ group, p.adjust.method = "holm")

dunn_res


## -------------------------
## 5C) 参考：pairwise Wilcoxon test
## Holm補正
## -------------------------
pairwise_wilcox_res <- pairwise.wilcox.test(
  x = df$value,
  g = df$group,
  p.adjust.method = "holm",
  exact = FALSE
)

pairwise_wilcox_res