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

library(tidyverse)
library(stringr)

# 12項目を定義
metrics_12 <- c(
  "Hip Angle range",
  "Knee Angle range",
  "Ankle Angle range",
  "Hip Angle Velocity Maximum",
  "Knee Angle Velocity Maximum",
  "Ankle Angle Velocity Maximum",
  "Hip Angle Velocity Minimum",
  "Knee Angle Velocity Minimum",
  "Ankle Angle Velocity Minimum",
  "Tip X Range",
  "Tip Y Range",
  "Tip Velocity XY Average"
)

device_7w$group    <- "Device"
ees_eeson_7w$group <- "EES_on"
ees_eesoff_7w$group <- "EES_off"

df_all <- bind_rows(
  device_7w,
  ees_eeson_7w,
  ees_eesoff_7w
)

# 1) ワイド→ロング（列名に埋まっている個体・群・左右を抜く）
df_long <- df_all %>%
  rename(metric = animal_name) %>%   # animal_name列は「metric」
  pivot_longer(
    cols = -c(metric, group),
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value)) %>%
  # raw_id 例: g2_device_ee_soff_7w_r
  separate(
    raw_id,
    into = c("animal", "cond1", "ee", "stim", "week", "side"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    animal = str_to_upper(animal),
    side = str_to_upper(side)
  )

# 12項目（データに実在する表記に合わせる）
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

# 1) 12項目だけ抽出
df_long_12 <- df_long %>%
  filter(metric %in% metrics_12_use)

# 念の為：ちゃんと12個そろっているかチェック
df_long_12 %>% distinct(metric) %>% arrange(metric)
df_long_12 %>% count(metric) %>% arrange(metric)

# 2) 左右平均（1 animal × group × metric で1つ）
df_meanLR <- df_long_12 %>%
  group_by(group, animal, metric) %>%
  summarise(val = mean(value, na.rm = TRUE), .groups = "drop")

# 3) PCA用ワイド化（1 animal=1行、12項目=12列）
df_pca <- df_meanLR %>%
  pivot_wider(names_from = metric, values_from = val) %>%
  drop_na()

# 4) PCA（スケールあり）
X <- df_pca %>% select(all_of(metrics_12_use))
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

# 5) スコア
pca_scores <- as_tibble(pca_res$x) %>%
  bind_cols(df_pca %>% select(group, animal))

# 6) 寄与率
summary(pca_res)

# 7) PC1-PC2 plot（95%楕円）
library(ggplot2)
ggplot(pca_scores, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(
    values = c(
      "EES_on"  = "#D62728",  # 赤
      "EES_off" = "#E377C2",  # ピンク
      "Device"  = "#4D4D4D"   # ダークグレー
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA (12 gait parameters)",
    x = "PC1 (60.4%)",
    y = "PC2 (15.4%)",
    color = "Group"
  )

library(vegan)

library(svglite)

ggsave(
  "~/Desktop/PCA_PC1_PC2_plot.pdf",
  plot = p,
  width = 6,
  height = 5
)

ggsave(
  "~/Desktop/PCA_PC1_PC2_plot.svg",
  plot = p,
  width = 6,
  height = 5
)

# PC1–PC2 行列
pc_mat <- pca_scores %>%
  select(PC1, PC2)

# ユークリッド距離
D <- dist(pc_mat, method = "euclidean")

set.seed(1)

# PERMANOVA（全体）
perm_res <- adonis2(
  D ~ group,
  data = pca_scores,
  permutations = 9999,
  by = "terms"
)

perm_res

# pairwise comparisons
groups <- combn(unique(pca_scores$group), 2, simplify = FALSE)

pairwise_perm <- purrr::map_dfr(groups, function(g) {
  
  df_sub <- pca_scores %>% filter(group %in% g)
  pc_sub <- df_sub %>% select(PC1, PC2)
  D_sub  <- dist(pc_sub)
  
  set.seed(1)
  res <- adonis2(
    D_sub ~ group,
    data = df_sub,
    permutations = 9999,
    by = "terms"
  )
  
  tibble(
    comparison = paste(g, collapse = " vs "),
    F = res$F[1],
    R2 = res$R2[1],
    p_raw = res$`Pr(>F)`[1]
  )
}) %>%
  mutate(
    p_BH = p.adjust(p_raw, method = "BH")
  ) %>%
  arrange(p_BH)

pairwise_perm

# PC1 データだけ抜き出す
pc1_df <- pca_scores %>%
  select(animal, group, PC1)

# 同一個体のみ抽出
pc1_on_off <- pc1_df %>%
  filter(group %in% c("EES_on", "EES_off")) %>%
  pivot_wider(names_from = group, values_from = PC1) %>%
  drop_na(EES_on, EES_off)

wilcox_on_off <- wilcox.test(
  pc1_on_off$EES_on,
  pc1_on_off$EES_off,
  paired = TRUE,
  exact = FALSE
)

wilcox_on_off

pc1_off_device <- pc1_df %>%
  filter(group %in% c("EES_off", "Device"))

wilcox_off_device <- wilcox.test(
  PC1 ~ group,
  data = pc1_off_device,
  exact = FALSE
)

wilcox_off_device

pc1_on_device <- pc1_df %>%
  filter(group %in% c("EES_on", "Device"))

wilcox_on_device <- wilcox.test(
  PC1 ~ group,
  data = pc1_on_device,
  exact = FALSE
)

wilcox_on_device

# p値まとめ
pvals <- c(
  EESon_vs_EESoff = wilcox_on_off$p.value,
  EESoff_vs_Device = wilcox_off_device$p.value,
  EESon_vs_Device = wilcox_on_device$p.value
)

pvals_adj <- p.adjust(pvals, method = "BH")

pvals
pvals_adj

# PCA と完全に同じ色
group_cols <- c(
  "EES_on"  = "#D62728",  # 赤
  "EES_off" = "#E377C2",  # ピンク
  "Device"  = "#4D4D4D"   # ダークグレー
)

# PC1用2群棒グラフ関数

library(tidyverse)

plot_pc1_bar <- function(pca_scores, g1, g2, paired = FALSE) {
  
  df <- pca_scores %>%
    filter(group %in% c(g1, g2)) %>%
    mutate(group = factor(group, levels = c(g1, g2))) %>%
    select(animal, group, PC1) %>%
    drop_na(PC1)
  
  # Mean ± SEM
  sum_df <- df %>%
    group_by(group) %>%
    summarise(
      mean = mean(PC1),
      sem  = sd(PC1) / sqrt(n()),
      .groups = "drop"
    )
  
  # 統計
  if (paired) {
    wide <- df %>%
      pivot_wider(names_from = group, values_from = PC1) %>%
      drop_na(!!sym(g1), !!sym(g2))
    
    p <- wilcox.test(
      wide[[g1]], wide[[g2]],
      paired = TRUE, exact = FALSE
    )$p.value
  } else {
    p <- wilcox.test(
      PC1 ~ group,
      data = df,
      exact = FALSE
    )$p.value
  }
  
  sig <- p < 0.05
  
  # 注釈位置
  y_max <- max(sum_df$mean + sum_df$sem)
  y_min <- min(sum_df$mean - sum_df$sem)
  y_bar <- y_max + 0.1 * (y_max - y_min + 1e-6)
  
  g <- ggplot(sum_df, aes(x = group, y = mean, fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    geom_col(width = 0.65, color = "black") +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width = 0.18, linewidth = 0.8
    ) +
    scale_fill_manual(values = group_cols) +   # ★ここが追加
    theme_classic(base_size = 15) +
    theme(legend.position = "none") +
    labs(
      title = paste0("PC1: ", g1, " vs ", g2),
      x = NULL,
      y = "PC1 score"
    )
  
  # ★ 有意なときだけ「＊」
  if (sig) {
    g <- g +
      geom_segment(
        aes(x = 1, xend = 2, y = y_bar, yend = y_bar),
        inherit.aes = FALSE, linewidth = 0.8
      ) +
      geom_text(
        aes(x = 1.5, y = y_bar, label = "*"),
        inherit.aes = FALSE, vjust = -0.4, size = 6
      )
  }
  
  g
}

plot_pc1_bar(pca_scores, "EES_on", "EES_off", paired = TRUE)

plot_pc1_bar(pca_scores, "EES_off", "Device", paired = FALSE)

plot_pc1_bar(pca_scores, "EES_on", "Device", paired = FALSE)