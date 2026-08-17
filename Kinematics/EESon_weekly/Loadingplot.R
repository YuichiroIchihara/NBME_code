############################################################
# PCA (12 gait parameters) + Loading plot (Arial, large text)
#  - Data: EES_EESon_weekly.xlsx
#  - Mean of L/R per animal
#  - PCA: prcomp(center=TRUE, scale.=TRUE)
#  - PC1 & PC2 are FLIPPED (scores + loadings)
#  - Loading plot: feature names shortened,
#                  ordered by |PC1| (TOP->BOTTOM)
############################################################


# ==========================================================
# 0. Packages
# ==========================================================

# 初回だけ必要
# install.packages(c("tidyverse", "readxl", "janitor", "stringr"))

library(tidyverse)
library(readxl)
library(janitor)
library(stringr)

# ==========================================================
# 1. Data read
# 添付の正常に動くスクリプトと同じpath
# ==========================================================

data_file <- "~/Desktop/Kinema_data/data/EES_EESon_weekly.xlsx"

# ファイルが存在するか確認
if (!file.exists(data_file)) {
  stop(
    paste0(
      "File not found: ",
      data_file
    )
  )
}

ees_eeson_weekly_raw <- read_excel(data_file)

# ==========================================================
# 2. clean_names()
# ==========================================================

ees_eeson_weekly <- ees_eeson_weekly_raw %>%
  clean_names()


# ==========================================================
# A) 「1行目がヘッダ」形式に対応
# ==========================================================

df <- ees_eeson_weekly

new_colnames <- as.character(df[1, ])

colnames(df) <- new_colnames

# 1行目（headerとして使用した行）を削除
df <- df[-1, ]

# 先頭列名を metric に統一
colnames(df)[1] <- "metric"

# ==========================================================
# B) long化して raw_id を分解
#
# raw_id例:
# g2_ees_eeson_7w_r
#
# animal / ees / ees_status / week / side
# ==========================================================

df_long <- df %>%
  pivot_longer(
    cols = -metric,
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value)
  ) %>%
  tidyr::separate(
    raw_id,
    into = c(
      "animal",
      "ees",
      "ees_status",
      "week",
      "side"
    ),
    sep = "_",
    remove = FALSE,
    fill = "right"
  ) %>%
  mutate(
    animal = str_to_upper(animal),
    side   = str_to_upper(side),
    week   = str_to_upper(week)
  )

# ==========================================================
# C) PCAに使用する12 metrics
# 最終的にこの名前に統一
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
# 表記ゆれを統一
# ==========================================================

df_long2 <- df_long %>%
  mutate(
    
    metric_std = case_when(
      
      # --------------------------
      # Angle range
      # --------------------------
      
      str_detect(
        metric,
        regex(
          "^Hip\\s*Angle\\s*range$",
          ignore_case = TRUE
        )
      ) ~ "Hip Angle range",
      
      str_detect(
        metric,
        regex(
          "^Knee\\s*Angle\\s*range$",
          ignore_case = TRUE
        )
      ) ~ "Knee Angle range",
      
      str_detect(
        metric,
        regex(
          "^Ankle\\s*Angle\\s*range$",
          ignore_case = TRUE
        )
      ) ~ "Ankle Angle range",
      
      
      # --------------------------
      # Velocity Maximum
      # --------------------------
      
      str_detect(
        metric,
        regex(
          "^Hip\\s*Angle\\s*Velocity\\s*Maximum$",
          ignore_case = TRUE
        )
      ) ~ "Hip Angle Velocity Maximum",
      
      str_detect(
        metric,
        regex(
          "^Knee\\s*Angle\\s*Velocity\\s*Maximum$",
          ignore_case = TRUE
        )
      ) ~ "Knee Angle Velocity Maximum",
      
      str_detect(
        metric,
        regex(
          "^Ankle\\s*Angle\\s*Velocity\\s*Maximum$",
          ignore_case = TRUE
        )
      ) ~ "Ankle Angle Velocity Maximum",
      
      
      # --------------------------
      # Velocity Minimum
      # → Minimum absolute に統一
      # --------------------------
      
      str_detect(
        metric,
        regex(
          "^Hip\\s*Angle\\s*Velocity\\s*Minimum(\\s*absolute)?$",
          ignore_case = TRUE
        )
      ) ~ "Hip Angle Velocity Minimum absolute",
      
      str_detect(
        metric,
        regex(
          "^Knee\\s*Angle\\s*Velocity\\s*Minimum(\\s*absolute)?$",
          ignore_case = TRUE
        )
      ) ~ "Knee Angle Velocity Minimum absolute",
      
      str_detect(
        metric,
        regex(
          "^Ankle\\s*Angle\\s*Velocity\\s*Minimum(\\s*absolute)?$",
          ignore_case = TRUE
        )
      ) ~ "Ankle Angle Velocity Minimum absolute",
      
      
      # --------------------------
      # Tip
      # --------------------------
      
      str_detect(
        metric,
        regex(
          "^Tip\\s*X\\s*Range$",
          ignore_case = TRUE
        )
      ) ~ "Tip X Range",
      
      str_detect(
        metric,
        regex(
          "^Tip\\s*Y\\s*Range$",
          ignore_case = TRUE
        )
      ) ~ "Tip Y Range",
      
      str_detect(
        metric,
        regex(
          "^Tip\\s*Velocity\\s*XY\\s*Average$",
          ignore_case = TRUE
        )
      ) ~ "Tip Velocity XY Average",
      
      TRUE ~ NA_character_
    ),
    
    # Velocity Minimum は絶対値化
    value2 = if_else(
      !is.na(metric_std) &
        str_detect(
          metric_std,
          "Velocity Minimum absolute"
        ),
      abs(value),
      value
    )
  )

# ==========================================================
# 12項目だけ残す
# ==========================================================

df_long_12 <- df_long2 %>%
  dplyr::filter(
    metric_std %in% metrics_12_use
  )

# ==========================================================
# CHECK
# 12種類すべて認識されているか
# ==========================================================

present_metrics <- df_long_12 %>%
  distinct(metric_std) %>%
  arrange(metric_std)

print(present_metrics)

cat(
  "\nNumber of metrics detected:",
  nrow(present_metrics),
  "\n"
)

# ==========================================================
# D) L/R 平均
#
# 1行 =
# animal × week × metric
# ==========================================================

df_meanLR <- df_long_12 %>%
  group_by(
    week,
    animal,
    metric = metric_std
  ) %>%
  summarise(
    val = mean(
      value2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ==========================================================
# E) PCA wide
#
# 1行 = animal × week
# 12列 = metrics
# ==========================================================

df_pca <- df_meanLR %>%
  pivot_wider(
    names_from = metric,
    values_from = val
  ) %>%
  drop_na(
    all_of(metrics_12_use)
  )


# データ確認
glimpse(df_pca)

# ==========================================================
# PCA入力行列
# ==========================================================

X <- df_pca %>%
  dplyr::select(
    all_of(metrics_12_use)
  )


# ==========================================================
# F) PCA
# ==========================================================

pca_res <- prcomp(
  X,
  center = TRUE,
  scale. = TRUE
)

# ==========================================================
# IMPORTANT
# PC1 & PC2 の向きを反転
#
# PCAの符号は任意なので、
# scoreとloadingを同時に反転する
# ==========================================================

pca_res$x[, c("PC1", "PC2")] <-
  -pca_res$x[, c("PC1", "PC2")]

pca_res$rotation[, c("PC1", "PC2")] <-
  -pca_res$rotation[, c("PC1", "PC2")]

# ==========================================================
# Explained variance
# ==========================================================

pve <- (
  pca_res$sdev^2 /
    sum(pca_res$sdev^2)
) * 100

cat(
  "PC1:",
  round(pve[1], 1),
  "%\n"
)

cat(
  "PC2:",
  round(pve[2], 1),
  "%\n"
)

# ==========================================================
# PCA scores
# ==========================================================

pca_scores <- as_tibble(
  pca_res$x
) %>%
  bind_cols(
    df_pca %>%
      dplyr::select(
        week,
        animal
      )
  )

# ==========================================================
# 10) PC1-PC2 scatter
# ==========================================================

p_pca <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = week
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  stat_ellipse(
    level = 0.95,
    linewidth = 1
  ) +
  
  scale_color_manual(
    values = c(
      "0W" = "#1f77b4",
      "3W" = "#2ca02c",
      "7W" = "#d62728"
    )
  ) +
  
  theme_classic(
    base_size = 14,
    base_family = "Arial"
  ) +
  
  labs(
    title = "PCA (12 gait parameters)",
    x = paste0(
      "PC1 (",
      round(pve[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pve[2], 1),
      "%)"
    ),
    color = "Week"
  )


print(p_pca)

# ==========================================================
# 11) Loading plot
# rename + |PC1|順 TOP → BOTTOM
# ==========================================================

rot <- as.data.frame(
  pca_res$rotation
) %>%
  tibble::rownames_to_column(
    "feature"
  )

# ==========================================================
# Feature表示名を短縮
# ==========================================================

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
# 12項目すべて使用
# ==========================================================

top_n <- 12


top_features <- rot2 %>%
  mutate(
    max_abs = pmax(
      abs(PC1),
      abs(PC2)
    )
  ) %>%
  arrange(
    desc(max_abs)
  ) %>%
  slice_head(
    n = top_n
  ) %>%
  pull(feature)

# ==========================================================
# Loading plot用long data
# ==========================================================

feature_order <- rot2 %>%
  dplyr::filter(
    feature %in% top_features
  ) %>%
  arrange(
    desc(abs(PC1))
  ) %>%
  pull(feature)


plot_df <- rot2 %>%
  dplyr::filter(
    feature %in% top_features
  ) %>%
  
  dplyr::select(
    feature,
    PC1,
    PC2
  ) %>%
  
  pivot_longer(
    cols = c(
      PC1,
      PC2
    ),
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
    
    feature = factor(
      feature,
      levels = feature_order
    )
  )

# ==========================================================
# Loading plot
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
  
  # |PC1| が大きいものを上から表示
  scale_y_discrete(
    limits = rev(
      levels(plot_df$feature)
    )
  ) +
  
  labs(
    title = "PCA Loading Plot (PC1 & PC2)",
    x = "Loading",
    y = "Feature",
    fill = "PC"
  ) +
  
  scale_fill_manual(
    values = c(
      "PC1" = "#8C510A",
      "PC2" = "#DFC27D"
    )
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
    
    legend.position = "right"
  )


print(p_loading)