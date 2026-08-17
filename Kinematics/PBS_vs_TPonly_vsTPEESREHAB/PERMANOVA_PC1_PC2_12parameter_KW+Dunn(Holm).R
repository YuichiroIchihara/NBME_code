install.packages(c(
  "tidyverse",  # dplyr, ggplot2 などまとめて入る便利セット
  "readxl",     # Excel読み込み用
  "janitor"     # 列名きれいにする用（あると便利）
))
library(tidyverse)
library(readxl)
library(janitor)

getwd()

ees_eeson_weekly_raw <- read_excel("~/Desktop/Kinema_data/data/EES_EESon_weekly.xlsx")
ees_eeson_7w_raw <- read_excel("~/Desktop/Kinema_data/data/EES_EESon_7W.xlsx")
ees_eesoff_7w_raw <- read_excel("~/Desktop/Kinema_data/data/EES_EESoff_7W.xlsx")
device_7w_raw <- read_excel("~/Desktop/Kinema_data/data/Device_7W.xlsx")
combination_7w_raw <- read_excel("~/Desktop/Kinema_data/data/Combination_7W.xlsx")
tp_7w_raw <- read_excel("~/Desktop/Kinema_data/data/TP_7W.xlsx")
pbs_7w_raw <- read_excel("~/Desktop/Kinema_data/data/PBS_7W.xlsx")

glimpse(ees_eeson_weekly_raw)
glimpse(ees_eeson_7w_raw)
glimpse(ees_eesoff_7w_raw)
glimpse(device_7w_raw)
glimpse(combination_7w_raw)
glimpse(tp_7w_raw)
glimpse(pbs_7w_raw)

ees_eeson_weekly <- ees_eeson_weekly_raw %>% clean_names()
ees_eeson_7w     <- ees_eeson_7w_raw     %>% clean_names()
ees_eesoff_7w    <- ees_eesoff_7w_raw    %>% clean_names()
device_7w        <- device_7w_raw        %>% clean_names()
combination_7w   <- combination_7w_raw   %>% clean_names()
tp_7w            <- tp_7w_raw            %>% clean_names()
pbs_7w           <- pbs_7w_raw           %>% clean_names()

# =========================
# ここから追加でOK
# =========================

# 0) もし read_excel をまだ skip=1 で読んでいない場合は読み直し（重要）
pbs_7w_raw <- read_excel("~/Desktop/Kinema_data/data/PBS_7W.xlsx", skip = 1)
tp_7w_raw  <- read_excel("~/Desktop/Kinema_data/data/TP_7W.xlsx",  skip = 1)
combination_7w_raw <- read_excel("~/Desktop/Kinema_data/data/Combination_7W.xlsx", skip = 1)

pbs_7w <- pbs_7w_raw %>% clean_names()
tp_7w  <- tp_7w_raw  %>% clean_names()
combination_7w <- combination_7w_raw %>% clean_names()

# =========================
# PBS / TP / Combination（7W）だけで 12項目PCA
# 前提：pbs_7w, tp_7w, combination_7w は clean_names() 済み
# =========================

library(tidyverse)
library(stringr)
library(janitor)

# 12項目（表記は「Excelの1列目（animal_name）に入っている文字列」に合わせる）
metrics_12_use <- c(
  "Hip Angle range",
  "Knee Angle range",
  "Ankle Angle range",
  "Hip Angle Velocity Maximum",
  "Knee Angle Velocity Maximum",
  "Ankle Angle Velocity Maximum",
  "Hip Angle Velocity Minimum absolute",
  "Knee Angle Velocity Minimum absolute",
  "Ankle Angle Velocity Minimum absolute",
  "Tip X Range",
  "Tip Y Range",
  "Tip Velocity XY Average"
)

# 1) group を付与して結合
pbs_7w$group         <- "PBS"
tp_7w$group          <- "TP"
combination_7w$group <- "Combination"

df_all <- bind_rows(pbs_7w, tp_7w, combination_7w)

# 2) ワイド→ロング（1列目 animal_name を metric として扱う）
df_long <- df_all %>%
  rename(metric = 1) %>%
  pivot_longer(
    cols = -c(metric, group),
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    # 先頭2トークンを animal にする：p3_1 → P3_1 （=別個体）
    animal = str_to_upper(str_extract(raw_id, "^[^_]+_[^_]+")),
    side   = str_to_upper(str_extract(raw_id, "[^_]+$"))
  )

# 3) 12項目だけ抽出
df_long_12 <- df_long %>%
  filter(metric %in% metrics_12_use)

# ---- チェック（12項目が揃っているか）----
present_metrics <- df_long_12 %>% distinct(metric) %>% pull(metric)
missing_metrics <- setdiff(metrics_12_use, present_metrics)
extra_metrics   <- setdiff(present_metrics, metrics_12_use)

cat("---- metric check ----\n")
cat("present:", length(present_metrics), "\n")
if (length(missing_metrics) > 0) cat("MISSING:\n", paste(missing_metrics, collapse = "\n"), "\n")
if (length(extra_metrics)   > 0) cat("EXTRA:\n",   paste(extra_metrics, collapse = "\n"), "\n")

# 4) 左右平均（1 animal × group × metric で1つ）
df_meanLR <- df_long_12 %>%
  group_by(group, animal, metric) %>%
  summarise(val = mean(value, na.rm = TRUE), .groups = "drop")

# 5) PCA用ワイド化（1 animal=1行、12項目=12列）
df_pca <- df_meanLR %>%
  pivot_wider(names_from = metric, values_from = val)

# NA がある個体は PCA から除外（必要なら imputations に変えてOK）
df_pca_complete <- df_pca %>%
  drop_na(all_of(metrics_12_use))

# 6) PCA（スケールあり）
X <- df_pca_complete %>% select(all_of(metrics_12_use))
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

# 7) スコア
pca_scores <- as_tibble(pca_res$x) %>%
  bind_cols(df_pca_complete %>% select(group, animal))

# 8) 寄与率（PC1/PC2% を自動でラベルに入れる）
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_pct <- round(100 * var_explained[1], 1)
pc2_pct <- round(100 * var_explained[2], 1)

print(summary(pca_res))

group_colors <- c(
  "PBS"         = "#0072B2",
  "TP"          = "#D55E00",
  "Combination" = "#2E8B57"
)

# =========================
# Plot（2枚目と同じデザインに寄せる）
# =========================

# 2枚目っぽいテーマを関数化（再利用可）
theme_pca_like2 <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      # グリッド（補助線）
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_line(linewidth = 0.25),
      
      # 枠線（パネルの外枠）
      panel.border = element_rect(linewidth = 0.8, fill = NA),
      
      # 軸線は消してOK（枠線があるので二重に見えやすい）
      axis.line = element_blank(),
      
      # 凡例など
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size * 0.9),
      plot.title   = element_text(hjust = 0)
    )
}

ggplot(pca_scores, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PCA (12 gait parameters)",
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)"),
    color = "Group"
  ) +
  theme_pca_like2(base_size = 14)

library(vegan)

# PC1–PC2 のみ抽出
pc12 <- pca_scores %>%
  select(PC1, PC2)

# ユークリッド距離
dist_pc12 <- dist(pc12, method = "euclidean")

perm_all <- adonis2(
  dist_pc12 ~ group,
  data = pca_scores,
  permutations = 999
)

print(perm_all)

# =========================
# PC1 の 3群比較：Kruskal-Wallis + Dunn
# =========================

# 必要パッケージ
install.packages(c("rstatix", "FSA"))  # 未インストールなら
library(rstatix)
library(FSA)
library(dplyr)
library(ggplot2)

# 0) PC1データ作成（pca_scores を使用）
pc1_df <- pca_scores %>%
  select(group, animal, PC1) %>%
  mutate(group = factor(group, levels = c("PBS", "TP", "Combination"))) %>%
  drop_na(PC1)

# 1) Kruskal-Wallis（全体差）
kw_res <- kruskal_test(pc1_df, PC1 ~ group)
kw_res

# 2) Dunn 検定（ペア比較）- BH補正（探索寄り）と Holm補正（堅め）を両方
dunn_bh <- dunnTest(PC1 ~ group, data = pc1_df, method = "bh")$res %>%
  as_tibble() %>%
  rename(p_adj = P.adj, z = Z) %>%
  mutate(p_raw = P.unadj, p_adj_method = "BH(FDR)") %>%
  select(comparison = Comparison, z, p_raw, p_adj, p_adj_method)

dunn_holm <- dunnTest(PC1 ~ group, data = pc1_df, method = "holm")$res %>%
  as_tibble() %>%
  rename(p_adj = P.adj, z = Z) %>%
  mutate(p_raw = P.unadj, p_adj_method = "Holm") %>%
  select(comparison = Comparison, z, p_raw, p_adj, p_adj_method)

dunn_res <- bind_rows(dunn_bh, dunn_holm) %>%
  arrange(p_adj_method, p_adj)

dunn_res

# 3) ついでに効果量（参考）：Kruskal-Wallis の epsilon^2
eff_kw <- kruskal_effsize(pc1_df, PC1 ~ group)  # rstatix
eff_kw

# 4) 図（棒グラフ or 箱ひげ：好みで）
# ---- 箱ひげ + 個体点（黒縁）----
ggplot(pc1_df, aes(x = group, y = PC1, fill = group)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color = group), width = 0.12, size = 2.6, stroke = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(title = "PC1 comparison (Kruskal-Wallis + Dunn)", x = NULL, y = "PC1 score") +
  theme_pca_like2(base_size = 14) +
  theme(legend.position = "none")


# ---- 棒（平均±SEM）+ 個体点（黒縁）----
pc1_sum <- pc1_df %>%
  group_by(group) %>%
  summarise(
    mean = mean(PC1),
    sem  = sd(PC1) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

ggplot() +
  geom_col(data = pc1_sum, aes(x = group, y = mean, fill = group),
           color = "black", width = 0.65, alpha = 0.85) +
  geom_errorbar(data = pc1_sum, aes(x = group, ymin = mean - sem, ymax = mean + sem),
                width = 0.18, linewidth = 0.8) +
  geom_jitter(data = pc1_df, aes(x = group, y = PC1, color = group),
              width = 0.12, size = 2.6, stroke = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(title = "PC1 comparison (mean ± SEM)", x = NULL, y = "PC1 score") +
  theme_pca_like2(base_size = 14) +
  theme(legend.position = "none")

# =========================
# PC1 図（mean±SEMの棒 + 個体点 + 0ライン）
# =========================

library(dplyr)
library(ggplot2)

# group の順序をここで固定（すでに pc1_df でfactor化済みなら不要）
pc1_df <- pc1_df %>%
  mutate(group = factor(group, levels = c("PBS", "TP", "Combination")))

# mean ± SEM
pc1_sum <- pc1_df %>%
  group_by(group) %>%
  summarise(
    mean = mean(PC1, na.rm = TRUE),
    sem  = sd(PC1, na.rm = TRUE) / sqrt(sum(!is.na(PC1))),
    n    = sum(!is.na(PC1)),
    .groups = "drop"
  )

# 色（PCAと同じにしたいなら、あなたの group_colors をそのまま使ってOK）
group_colors <- c(
  "PBS"         = "#0072B2",
  "TP"          = "#D55E00",
  "Combination" = "#2E8B57"
)

p_pc1_bar <- ggplot() +
  # 0ライン
  geom_hline(yintercept = 0, linewidth = 0.9) +
  # バー（黒枠）
  geom_col(
    data = pc1_sum,
    aes(x = group, y = mean, fill = group),
    color = "black",
    width = 0.62
  ) +
  # SEM
  geom_errorbar(
    data = pc1_sum,
    aes(x = group, ymin = mean - sem, ymax = mean + sem),
    width = 0.18,
    linewidth = 0.9
  ) +
  # 個体点（黒縁・塗りは群色）
  geom_jitter(
    data = pc1_df,
    aes(x = group, y = PC1, fill = group),
    shape = 21,
    color = "black",
    size = 3,
    stroke = 0.7,
    width = 0.08
  ) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "PC1 score") +
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_line(linewidth = 1.1),
    axis.line.y = element_line(linewidth = 1.1),
    axis.ticks  = element_line(linewidth = 1.1),
    legend.position = "none"
  )

p_pc1_bar

# =========================
# PC2 の 3群比較：Kruskal-Wallis + Dunn
# =========================

install.packages(c("rstatix", "FSA"))
library(rstatix)
library(FSA)
library(dplyr)
library(ggplot2)
library(tidyr)

# 0) PC2データ作成
pc2_df <- pca_scores %>%
  select(group, animal, PC2) %>%
  mutate(group = factor(group, levels = c("PBS", "TP", "Combination"))) %>%
  drop_na(PC2)

# 1) Kruskal-Wallis（全体差）
kw_res_pc2 <- kruskal_test(pc2_df, PC2 ~ group)
kw_res_pc2

# =========================
# PC2 図（0 baseline の棒 + 個体点）
# =========================

# mean ± SEM
pc2_sum <- pc2_df %>%
  group_by(group) %>%
  summarise(
    mean = mean(PC2, na.rm = TRUE),
    sem  = sd(PC2, na.rm = TRUE) / sqrt(sum(!is.na(PC2))),
    n    = sum(!is.na(PC2)),
    .groups = "drop"
  )

p_pc2_bar <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.9) +
  geom_col(
    data = pc2_sum,
    aes(x = group, y = mean, fill = group),
    color = "black",
    width = 0.62
  ) +
  geom_errorbar(
    data = pc2_sum,
    aes(x = group, ymin = mean - sem, ymax = mean + sem),
    width = 0.18,
    linewidth = 0.9
  ) +
  geom_jitter(
    data = pc2_df,
    aes(x = group, y = PC2, fill = group),
    shape = 21,
    color = "black",
    size = 3,
    stroke = 0.7,
    width = 0.08
  ) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "PC2 score") +
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_line(linewidth = 1.1),
    axis.line.y = element_line(linewidth = 1.1),
    axis.ticks  = element_line(linewidth = 1.1),
    legend.position = "none"
  )

p_pc2_bar

# =========================
# PERMANOVA (PC1–PC2)
# =========================

install.packages("vegan")
library(vegan)
library(dplyr)

# PCAスコア（PC1, PC2）だけ抽出
permanova_df <- pca_scores %>%
  select(group, animal, PC1, PC2) %>%
  mutate(group = factor(group, levels = c("PBS", "TP", "Combination"))) %>%
  drop_na(PC1, PC2)

# 距離行列（ユークリッド距離）
dist_pc12 <- dist(permanova_df[, c("PC1", "PC2")], method = "euclidean")

# PERMANOVA（全体）
set.seed(123)  # 再現性
adonis_pc12 <- adonis2(
  dist_pc12 ~ group,
  data = permanova_df,
  permutations = 999
)

adonis_pc12

# =========================
# 12項目それぞれの箱ひげ図 + 統計
# Kruskal-Wallis + Dunn's test with Holm correction
# =========================

# 必要パッケージ
if (!requireNamespace("rstatix", quietly = TRUE)) {
  install.packages("rstatix")
}

library(tidyverse)
library(rstatix)
library(ggplot2)

# group順、metric順を固定
metric_df <- df_meanLR %>%
  mutate(
    group = factor(group, levels = c("PBS", "TP", "Combination")),
    metric = factor(metric, levels = metrics_12_use)
  ) %>%
  drop_na(val)

# 念のため確認
metric_df %>%
  count(metric, group)

# =========================
# 1) 各metricの記述統計
# =========================

metric_summary <- metric_df %>%
  group_by(metric, group) %>%
  summarise(
    n = n(),
    mean = mean(val, na.rm = TRUE),
    sd = sd(val, na.rm = TRUE),
    sem = sd / sqrt(n),
    median = median(val, na.rm = TRUE),
    min = min(val, na.rm = TRUE),
    max = max(val, na.rm = TRUE),
    .groups = "drop"
  )

print(metric_summary)

# =========================
# 2) 各metricごとのKruskal-Wallis検定
# =========================

kw_metric <- metric_df %>%
  group_by(metric) %>%
  kruskal_test(val ~ group) %>%
  arrange(p)

print(kw_metric)

# =========================
# 3) Dunn検定 Holm補正
# =========================

dunn_metric <- metric_df %>%
  group_by(metric) %>%
  dunn_test(val ~ group, p.adjust.method = "holm") %>%
  arrange(metric, p.adj)

print(dunn_metric)