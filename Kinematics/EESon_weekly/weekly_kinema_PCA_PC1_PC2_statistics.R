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

glimpse(ees_eeson_weekly)

library(tidyverse)

df <- ees_eeson_weekly

# --- 1行目を列名として設定 ---
new_colnames <- as.character(df[1, ])
colnames(df) <- new_colnames

# 1行目（ヘッダだった場所）を削除
df <- df[-1, ]

# 列名を整える
colnames(df)[1] <- "metric"

# --- ワイド → ロングに再変換 ---
df_long <- df %>%
  pivot_longer(
    cols = -metric,
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value))

# --- animal / week / side を抽出 ---
df_long <- df_long %>%
  separate(raw_id, into = c("animal", "ees", "ees_status", "week", "side"), sep = "_")

pca_long <- df_long %>%
  filter(
    str_detect(metric, regex(
      "Hip Angle Range|
       Knee Angle Range|
       Ankle Angle Range|
       Hip Angle Velocity Maximum|
       Knee Angle Velocity Maximum|
       Ankle Angle Velocity Maximum|
       Hip Angle Velocity Minimum|
       Knee Angle Velocity Minimum|
       Ankle Angle Velocity Minimum|
       Tip X Range|
       Tip Y Range|
       Tip Velocity XY Average",
      ignore_case = TRUE
    ))
  ) %>%
  mutate(
    value2 = case_when(
      str_detect(metric, regex("Velocity Minimum", ignore_case = TRUE)) ~ abs(value),
      TRUE ~ value
    )
  ) %>%
  group_by(animal, week, metric) %>%
  summarise(val = mean(value2, na.rm = TRUE), .groups = "drop")

pca_df <- pca_long %>%
  pivot_wider(
    names_from  = metric,
    values_from = val
  ) %>%
  drop_na()

glimpse(pca_df)

# df_long がある前提

# range が存在するか
df_long %>%
  distinct(metric) %>%
  filter(str_detect(metric, regex("Angle\\s*range", ignore_case = TRUE))) %>%
  print(n = 200)

# Tip X/Y Range が存在するか
df_long %>%
  distinct(metric) %>%
  filter(str_detect(metric, regex("^Tip\\s*X\\s*Range|^Tip\\s*Y\\s*Range", ignore_case = TRUE))) %>%
  print(n = 200)

# 12項目に正規化した名前を作る
df_long2 <- df_long %>%
  mutate(
    metric_std = case_when(
      # --- Angle Range (3) ---
      str_detect(metric, regex("^Hip\\s*Angle\\s*range$",  TRUE)) ~ "Hip Angle Range",
      str_detect(metric, regex("^Knee\\s*Angle\\s*range$", TRUE)) ~ "Knee Angle Range",
      str_detect(metric, regex("^Ankle\\s*Angle\\s*range$",TRUE)) ~ "Ankle Angle Range",
      
      str_detect(metric, regex("^Hip\\s*Angle\\s*Range$",  TRUE)) ~ "Hip Angle Range",
      str_detect(metric, regex("^Knee\\s*Angle\\s*Range$", TRUE)) ~ "Knee Angle Range",
      str_detect(metric, regex("^Ankle\\s*Angle\\s*Range$",TRUE)) ~ "Ankle Angle Range",
      
      # --- Velocity Max (3) ---
      str_detect(metric, regex("^Hip\\s*Angle\\s*Velocity\\s*Maximum$",  TRUE)) ~ "Hip Angle Velocity Maximum",
      str_detect(metric, regex("^Knee\\s*Angle\\s*Velocity\\s*Maximum$", TRUE)) ~ "Knee Angle Velocity Maximum",
      str_detect(metric, regex("^Ankle\\s*Angle\\s*Velocity\\s*Maximum$",TRUE)) ~ "Ankle Angle Velocity Maximum",
      
      # --- Velocity Min (3) ※ absolute化する ---
      str_detect(metric, regex("^Hip\\s*Angle\\s*Velocity\\s*Minimum(\\s*absolute)?$",  TRUE)) ~ "Hip Angle Velocity Minimum",
      str_detect(metric, regex("^Knee\\s*Angle\\s*Velocity\\s*Minimum(\\s*absolute)?$", TRUE)) ~ "Knee Angle Velocity Minimum",
      str_detect(metric, regex("^Ankle\\s*Angle\\s*Velocity\\s*Minimum(\\s*absolute)?$",TRUE)) ~ "Ankle Angle Velocity Minimum",
      
      # --- Tip (3) ---
      str_detect(metric, regex("^Tip\\s*X\\s*Range$", TRUE)) ~ "Tip X Range",
      str_detect(metric, regex("^Tip\\s*Y\\s*Range$", TRUE)) ~ "Tip Y Range",
      str_detect(metric, regex("^Tip\\s*Velocity\\s*XY\\s*Average$", TRUE)) ~ "Tip Velocity XY Average",
      
      TRUE ~ NA_character_
    )
  )

pca_long <- df_long2 %>%
  filter(!is.na(metric_std)) %>%
  mutate(
    value2 = if_else(
      str_detect(metric_std, "Velocity Minimum"),
      abs(value),
      value
    )
  ) %>%
  group_by(animal, week, metric_std) %>%
  summarise(val = mean(value2, na.rm = TRUE), .groups = "drop")

# 何個の指標が揃ったか確認（ここが12になってほしい）
pca_long %>% distinct(metric_std) %>% arrange(metric_std) %>% print(n = 50)

pca_df <- pca_long %>%
  pivot_wider(names_from = metric_std, values_from = val)

glimpse(pca_df)

library(tidyverse)
library(ggplot2)
library(vegan)

# 解析行列（12項目のみ）
X <- pca_df %>% select(-animal, -week)

# PCA（z-score 標準化）
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

# 寄与率（%）
pve <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
pve[1:5]  # まずはPC1-5を確認

scores <- as.data.frame(pca_res$x) %>%
  mutate(animal = pca_df$animal,
         week   = pca_df$week)

ggplot(scores, aes(PC1, PC2, color = week)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  labs(
    title = "PCA (12 gait parameters)",
    x = paste0("PC1 (", round(pve[1], 1), "%)"),
    y = paste0("PC2 (", round(pve[2], 1), "%)")
  ) +
  scale_color_manual(
    values = c("0W"="#1f77b4","3W"="#2ca02c","7W"="#d62728")
  ) +
  theme_bw(base_size = 14)

# =========================================
# PERMANOVA on PC1–PC2 distances (repeated measures)
#  - distance: Euclidean on (PC1, PC2)
#  - factor: week
#  - permutation restricted within animal (strata)
# =========================================
library(vegan)

# 解析に使うPC（PC1, PC2）
pc_mat <- scores %>%
  select(PC1, PC2)

# 距離行列（PC1–PC2のユークリッド距離）
D <- dist(pc_mat, method = "euclidean")

# week の順序を固定（任意）
scores <- scores %>%
  mutate(week = factor(week, levels = c("0W", "3W", "7W"), ordered = TRUE))

set.seed(1)

# PERMANOVA（個体内で置換を制限）
perm_res <- adonis2(
  D ~ week,
  data = scores,
  permutations = 9999,
  strata = scores$animal,   # ★ここが反復測定の肝
  by = "terms"
)

perm_res

# PC1データの準備
pc1_df <- scores %>%
  select(animal, week, PC1) %>%
  mutate(
    week = factor(week, levels = c("0W", "3W", "7W"), ordered = TRUE)
  )

# 全週そろっている個体のみ
pc1_wide <- pc1_df %>%
  pivot_wider(names_from = week, values_from = PC1) %>%
  drop_na(`0W`, `3W`, `7W`)

# Friedman検定
friedman.test(
  as.matrix(pc1_wide[, c("0W", "3W", "7W")])
)

# =========================================
# PC1 spaghetti plot
# =========================================
library(tidyverse)

pc1_long <- scores %>%
  select(animal, week, PC1) %>%
  mutate(
    week = factor(week, levels = c("0W", "3W", "7W"), ordered = TRUE)
  )

# 全週そろっている個体のみ（Friedmanと整合）
pc1_long <- pc1_long %>%
  group_by(animal) %>%
  filter(all(c("0W","3W","7W") %in% week)) %>%
  ungroup()

# Mean ± SEM
pc1_summary <- pc1_long %>%
  group_by(week) %>%
  summarise(
    mean = mean(PC1),
    sem  = sd(PC1) / sqrt(n()),
    .groups = "drop"
  )

ggplot(pc1_long, aes(x = week, y = PC1)) +
  geom_line(aes(group = animal), color = "grey70", linewidth = 0.6) +
  geom_point(aes(color = week), size = 3) +
  
  geom_line(
    data = pc1_summary,
    mapping = aes(x = week, y = mean, group = 1),
    inherit.aes = FALSE,
    color = "orange",
    linewidth = 1.2
  ) +
  geom_point(
    data = pc1_summary,
    mapping = aes(x = week, y = mean),
    inherit.aes = FALSE,
    color = "orange",
    size = 4
  ) +
  geom_errorbar(
    data = pc1_summary,
    mapping = aes(x = week, ymin = mean - sem, ymax = mean + sem),
    inherit.aes = FALSE,
    width = 0.15,
    color = "orange"
  ) +
  
  scale_color_manual(values = c("0W"="#1f77b4","3W"="#2ca02c","7W"="#d62728")) +
  labs(
    title = "PC1 score over time",
    y = "PC1 score",
    x = "Week"
  ) +
  theme_classic(base_size = 13)

# =========================================
# PC2 data (long + summary)
# =========================================
pc2_long <- scores %>%
  select(animal, week, PC2) %>%
  mutate(
    week = factor(week, levels = c("0W", "3W", "7W"), ordered = TRUE)
  )

# 全週そろっている個体のみ
pc2_long <- pc2_long %>%
  group_by(animal) %>%
  filter(all(c("0W","3W","7W") %in% week)) %>%
  ungroup()

pc2_summary <- pc2_long %>%
  group_by(week) %>%
  summarise(
    mean = mean(PC2),
    sem  = sd(PC2) / sqrt(n()),
    .groups = "drop"
  )

ggplot(pc2_long, aes(x = week, y = PC2)) +
  geom_line(aes(group = animal), color = "grey70", linewidth = 0.6) +
  geom_point(aes(color = week), size = 3) +
  
  geom_line(
    data = pc2_summary,
    mapping = aes(x = week, y = mean, group = 1),
    inherit.aes = FALSE,
    color = "orange",
    linewidth = 1.2
  ) +
  geom_point(
    data = pc2_summary,
    mapping = aes(x = week, y = mean),
    inherit.aes = FALSE,
    color = "orange",
    size = 4
  ) +
  geom_errorbar(
    data = pc2_summary,
    mapping = aes(x = week, ymin = mean - sem, ymax = mean + sem),
    inherit.aes = FALSE,
    width = 0.15,
    color = "orange"
  ) +
  
  scale_color_manual(values = c("0W"="#1f77b4","3W"="#2ca02c","7W"="#d62728")) +
  labs(
    title = "PC2 score over time",
    y = "PC2 score",
    x = "Week"
  ) +
  theme_classic(base_size = 13)

pc2_wide <- pc2_long %>%
  pivot_wider(names_from = week, values_from = PC2) %>%
  drop_na(`0W`, `3W`, `7W`)

friedman.test(
  as.matrix(pc2_wide[, c("0W", "3W", "7W")])
)