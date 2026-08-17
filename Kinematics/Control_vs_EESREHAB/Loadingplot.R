############################################################
# PCA (12 gait parameters) + Loading plot (Arial, large text)
#  - Groups: Device / EES_on / EES_off
#  - Mean of L/R per animal
#  - PCA: prcomp(center=TRUE, scale.=TRUE)
#  - Loading plot: feature names shortened,
#                  ordered by |PC1| (top -> bottom)
############################################################


# ==========================================================
# 0) Packages
# ==========================================================

# 初回だけ必要。すでに入っていれば実行しない
# install.packages(c("tidyverse", "readxl", "janitor", "stringr"))

library(tidyverse)
library(readxl)
library(janitor)
library(stringr)

# ==========================================================
# 1) Data read
# 正常に回っているスクリプトと同じデータ場所
# ==========================================================

ees_eeson_weekly_raw <- read_excel(
  "~/Desktop/Kinema_data/data/EES_EESon_weekly.xlsx"
)

ees_eeson_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/EES_EESon_7W.xlsx"
)

ees_eesoff_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/EES_EESoff_7W.xlsx"
)

device_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/Device_7W.xlsx"
)

combination_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/Combination_7W.xlsx"
)

tp_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/TP_7W.xlsx"
)

pbs_7w_raw <- read_excel(
  "~/Desktop/Kinema_data/data/PBS_7W.xlsx"
)


# ==========================================================
# Optional check
# ==========================================================

glimpse(ees_eeson_weekly_raw)
glimpse(ees_eeson_7w_raw)
glimpse(ees_eesoff_7w_raw)
glimpse(device_7w_raw)
glimpse(combination_7w_raw)
glimpse(tp_7w_raw)
glimpse(pbs_7w_raw)


# ==========================================================
# 2) clean_names()
# ==========================================================

ees_eeson_weekly <- ees_eeson_weekly_raw %>%
  clean_names()

ees_eeson_7w <- ees_eeson_7w_raw %>%
  clean_names()

ees_eesoff_7w <- ees_eesoff_7w_raw %>%
  clean_names()

device_7w <- device_7w_raw %>%
  clean_names()

combination_7w <- combination_7w_raw %>%
  clean_names()

tp_7w <- tp_7w_raw %>%
  clean_names()

pbs_7w <- pbs_7w_raw %>%
  clean_names()

# ==========================================================
# 3) PCAに使用する12 metrics
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
# 4) Group labels
# ==========================================================

device_7w$group <- "Device"
ees_eeson_7w$group <- "EES_on"
ees_eesoff_7w$group <- "EES_off"


# 3群を結合
df_all <- bind_rows(
  device_7w,
  ees_eeson_7w,
  ees_eesoff_7w
)

# ==========================================================
# 5) Wide -> Long
#
# animal_name列にはmetric名が入っているので
# metricにrename
# ==========================================================

df_long <- df_all %>%
  
  dplyr::rename(
    metric = animal_name
  ) %>%
  
  pivot_longer(
    cols = -c(metric, group),
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  
  mutate(
    value = as.numeric(value)
  ) %>%
  
  # raw_id example:
  # g2_device_ee_soff_7w_r
  tidyr::separate(
    raw_id,
    into = c(
      "animal",
      "cond1",
      "ee",
      "stim",
      "week",
      "side"
    ),
    sep = "_",
    remove = FALSE
  ) %>%
  
  mutate(
    animal = str_to_upper(animal),
    side = str_to_upper(side)
  )

# ==========================================================
# 6) 12 metricsだけ抽出
# ==========================================================

df_long_12 <- df_long %>%
  dplyr::filter(
    metric %in% metrics_12_use
  )


# ==========================================================
# Sanity check
# 12項目揃っているか確認
# ==========================================================

df_long_12 %>%
  distinct(metric) %>%
  arrange(metric) %>%
  print(n = 20)

df_long_12 %>%
  count(metric) %>%
  arrange(metric) %>%
  print(n = 20)

# ==========================================================
# 7) L/R平均
#
# 1 animal × group × metric = 1 value
# ==========================================================

df_meanLR <- df_long_12 %>%
  
  group_by(
    group,
    animal,
    metric
  ) %>%
  
  summarise(
    val = mean(
      value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ==========================================================
# 8) PCA用wide table
#
# 1 animal = 1 row
# 12 metrics = 12 columns
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
# 9) PCA
# ==========================================================

X <- df_pca %>%
  dplyr::select(
    all_of(metrics_12_use)
  )


pca_res <- prcomp(
  X,
  center = TRUE,
  scale. = TRUE
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
        group,
        animal
      )
  )

# ==========================================================
# Explained variance
# ==========================================================

summary(pca_res)


# PC1 / PC2 の寄与率を自動取得
pve <- (
  pca_res$sdev^2 /
    sum(pca_res$sdev^2)
) * 100

cat(
  "PC1 =",
  round(pve[1], 1),
  "%\n"
)

cat(
  "PC2 =",
  round(pve[2], 1),
  "%\n"
)

# ==========================================================
# 10) PC1-PC2 scatter plot
# ==========================================================

p_pca <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = group
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
      "EES_on"  = "#D62728",
      "EES_off" = "#E377C2",
      "Device"  = "#4D4D4D"
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
    
    color = "Group"
  )


print(p_pca)

############################################################
# 11) Loading plot
#     rename + order by |PC1|
#     TOP -> BOTTOM
############################################################


# ==========================================================
# Loading table
# ==========================================================

rot <- as.data.frame(
  pca_res$rotation
) %>%
  
  tibble::rownames_to_column(
    "feature"
  )

# ==========================================================
# Feature名を短縮
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
# 12項目すべて表示
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
# |PC1|順
# ==========================================================

feature_order <- rot2 %>%
  
  dplyr::filter(
    feature %in% top_features
  ) %>%
  
  arrange(
    desc(abs(PC1))
  ) %>%
  
  pull(feature)

# ==========================================================
# Loading plot用データ
# ==========================================================

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
  
  # |PC1|が大きいものを上から表示
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