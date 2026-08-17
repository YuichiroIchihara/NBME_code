# ==========================================================
# PBS / TP / Combination (7W)
# 12 gait parameters PCA + loading plot
# ==========================================================

library(tidyverse)
library(readxl)
library(janitor)
library(stringr)

# ==========================================================
# 1. データ読み込み
# 添付スクリプトと同じ場所・ファイル名を使用
# ==========================================================

pbs_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/PBS_7W.xlsx",
  skip = 1
)

tp_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/TP_7W.xlsx",
  skip = 1
)

combination_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/Combination_7W.xlsx",
  skip = 1
)

# 列名を整理
pbs_7w <- pbs_7w_raw %>% clean_names()
tp_7w <- tp_7w_raw %>% clean_names()
combination_7w <- combination_7w_raw %>% clean_names()

# ==========================================================
# 2. PCAに使用する12項目
# ==========================================================

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

# ==========================================================
# 3. groupを付けて3群を結合
# ==========================================================

pbs_7w$group <- "PBS"
tp_7w$group <- "TP"
combination_7w$group <- "Combination"

df_all <- bind_rows(
  pbs_7w,
  tp_7w,
  combination_7w
)

# ==========================================================
# 4. wide → long
# ==========================================================

df_long <- df_all %>%
  dplyr::rename(metric = 1) %>%
  pivot_longer(
    cols = -c(metric, group),
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    
    # 例：p3_1_r → animal = P3_1
    animal = str_to_upper(
      str_extract(raw_id, "^[^_]+_[^_]+")
    ),
    
    # 最後の部分をL/Rとして抽出
    side = str_to_upper(
      str_extract(raw_id, "[^_]+$")
    )
  )

# ==========================================================
# 5. 12項目だけ抽出
# ==========================================================

df_long_12 <- df_long %>%
  dplyr::filter(metric %in% metrics_12_use)

# ==========================================================
# 6. 12項目が正しく存在するか確認
# ==========================================================

present_metrics <- df_long_12 %>%
  distinct(metric) %>%
  pull(metric)

missing_metrics <- setdiff(
  metrics_12_use,
  present_metrics
)

extra_metrics <- setdiff(
  present_metrics,
  metrics_12_use
)

cat("---- metric check ----\n")
cat("present:", length(present_metrics), "\n")

if (length(missing_metrics) > 0) {
  cat(
    "MISSING:\n",
    paste(missing_metrics, collapse = "\n"),
    "\n"
  )
}

if (length(extra_metrics) > 0) {
  cat(
    "EXTRA:\n",
    paste(extra_metrics, collapse = "\n"),
    "\n"
  )
}

# ==========================================================
# 7. 左右平均
# 1 animal × group × metric = 1 value
# ==========================================================

df_meanLR <- df_long_12 %>%
  group_by(group, animal, metric) %>%
  summarise(
    val = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================================
# 8. PCA用にwide化
# 1 animal = 1行
# ==========================================================

df_pca <- df_meanLR %>%
  pivot_wider(
    names_from = metric,
    values_from = val
  )


# PCAに必要な12項目にNAがある個体は除外
df_pca_complete <- df_pca %>%
  drop_na(all_of(metrics_12_use))

# ==========================================================
# 9. PCA
# ==========================================================

X <- df_pca_complete %>%
  dplyr::select(all_of(metrics_12_use))

pca_res <- prcomp(
  X,
  center = TRUE,
  scale. = TRUE
)

# ==========================================================
# 10. PCA score
# ==========================================================

pca_scores <- as_tibble(pca_res$x) %>%
  bind_cols(
    df_pca_complete %>%
      dplyr::select(group, animal)
  )

# ==========================================================
# 11. PC1 / PC2 寄与率
# ==========================================================

var_explained <- (pca_res$sdev^2) /
  sum(pca_res$sdev^2)

pc1_pct <- round(
  100 * var_explained[1],
  1
)

pc2_pct <- round(
  100 * var_explained[2],
  1
)

print(summary(pca_res))

# ==========================================================
# 12. PCA plot
# ==========================================================

group_colors <- c(
  "PBS" = "#0072B2",
  "TP" = "#D55E00",
  "Combination" = "#2E8B57"
)

p_pca <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = group
  )
) +
  geom_point(size = 3) +
  stat_ellipse(
    level = 0.95,
    linewidth = 1
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  theme_classic(
    base_size = 14,
    base_family = "Arial"
  ) +
  labs(
    title = "PCA (12 gait parameters) : PBS vs TP vs Combination (7W)",
    x = paste0(
      "PC1 (",
      pc1_pct,
      "%)"
    ),
    y = paste0(
      "PC2 (",
      pc2_pct,
      "%)"
    ),
    color = "Group"
  )
print(p_pca)

# ==========================================================
# 13. Loading plot
# ==========================================================

rot <- as.data.frame(
  pca_res$rotation
) %>%
  tibble::rownames_to_column(
    "feature"
  ) %>%
  dplyr::filter(
    feature %in% metrics_12_use
  )

# 表示名を短くする
feature_map <- c(
  "Ankle Angle Velocity Maximum" =
    "Ankle Vel Max(Ext)",
  
  "Ankle Angle Velocity Minimum absolute" =
    "Ankle Vel Max(Flex)",
  
  "Hip Angle range" =
    "Hip Range",
  
  "Knee Angle Velocity Maximum" =
    "Knee Vel Max(Ext)",
  
  "Knee Angle range" =
    "Knee Range",
  
  "Tip Velocity XY Average" =
    "Tip Vel XY Avg",
  
  "Hip Angle Velocity Minimum absolute" =
    "Hip Vel Max(Flex)",
  
  "Hip Angle Velocity Maximum" =
    "Hip Vel Max(Ext)",
  
  "Knee Angle Velocity Minimum absolute" =
    "Knee Vel Max(Flex)",
  
  "Ankle Angle range" =
    "Ankle Range",
  
  "Tip X Range" =
    "Tip X Range",
  
  "Tip Y Range" =
    "Tip Y Range"
)
rot2 <- rot %>%
  mutate(
    feature_raw = feature,
    feature = dplyr::recode(
      feature,
      !!!feature_map,
      .default = feature
    )
  )

# ==========================================================
# 14. Loading plot用データ
# PC1 / PC2をlong形式に
# ==========================================================

plot_df <- rot2 %>%
  dplyr::select(
    feature,
    PC1,
    PC2
  ) %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "loading"
  ) %>%
  mutate(
    PC = factor(
      PC,
      levels = c(
        "PC1",
        "PC2"
      )
    ),
    
    # |PC1 loading| が大きい順
    feature = factor(
      feature,
      levels = rot2 %>%
        arrange(
          desc(abs(PC1))
        ) %>%
        pull(feature)
    )
  )

# ==========================================================
# 15. Loading plot
# ==========================================================

p_loading <- ggplot(
  plot_df,
  aes(
    x = loading,
    y = feature,
    fill = PC
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.65,
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.7
  ) +
  scale_y_discrete(
    limits = rev(
      levels(plot_df$feature)
    )
  ) +
  scale_fill_manual(
    values = c(
      "PC1" = "#8C510A",
      "PC2" = "#DFC27D"
    )
  ) +
  labs(
    title = "PCA Loading Plot (PC1 & PC2)",
    x = "Loading",
    y = "Feature",
    fill = "PC"
  ) +
  theme_bw(
    base_size = 18,
    base_family = "Arial"
  ) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 20
    ),
    axis.title.y = element_text(
      size = 20
    ),
    axis.text.x = element_text(
      size = 16
    ),
    axis.text.y = element_text(
      size = 18
    ),
    legend.title = element_text(
      size = 18
    ),
    legend.text = element_text(
      size = 16
    ),
    legend.position = "right",
    panel.grid.minor =
      element_blank()
  )
print(p_loading)