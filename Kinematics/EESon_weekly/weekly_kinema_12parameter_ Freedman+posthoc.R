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

# Ankle Angle range
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

# --- Ankle Angle Range 行だけ抽出 ---
ankle_range <- df_long %>%
  filter(str_detect(metric, regex("Ankle Angle range", ignore_case = TRUE)))

# 結果を確認
glimpse(ankle_range)

# 左右の平均化
ankle_range_mean <- ankle_range %>%
  group_by(animal, week) %>%
  summarise(
    ankle_angle_range_mean = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# 週の順序を指定(0W→3W→7W)
ankle_range_mean <- ankle_range_mean %>%
  mutate(week = factor(week, levels = c("0W", "3W", "7W")))

library(tidyverse)

# =========================================================
# Friedman (overall) + posthoc (paired) : Ankle Angle range
# =========================================================

# 1) 3時点（0W/3W/7W）が揃っている個体だけに絞る
ankle_complete <- ankle_range_mean %>%
  filter(week %in% c("0W","3W","7W")) %>%
  group_by(animal) %>%
  filter(sum(!is.na(ankle_angle_range_mean)) == 3) %>%
  ungroup()

# 2) ワイドにしてFriedman
ankle_wide <- ankle_complete %>%
  select(animal, week, ankle_angle_range_mean) %>%
  pivot_wider(names_from = week, values_from = ankle_angle_range_mean) %>%
  drop_na()

# --- Friedman test ---
friedman_ankle <- friedman.test(
  as.matrix(ankle_wide %>% select(`0W`, `3W`, `7W`))
)
friedman_ankle

# =========================================================
# 12 gait parameters:
# 左右平均 → 0W/3W/7W complete animal抽出 → Friedman test
# optional: paired posthoc Wilcoxon
# =========================================================

library(tidyverse)

# --- ワイド → ロングに再変換 ---
df_long <- df %>%
  pivot_longer(
    cols = -metric,
    names_to = "raw_id",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    
    # raw_idのどこかにある 0W / 3W / 7W を抽出
    week = str_extract(raw_id, regex("0W|3W|7W", ignore_case = TRUE)),
    week = str_to_upper(week),
    
    # sideを抽出
    side = str_extract(raw_id, regex("left|right|l|r$", ignore_case = TRUE)),
    side = str_to_upper(side),
    
    # animal名を作成
    # まず week と side 以降を除去
    animal = raw_id %>%
      str_remove(regex("_(0W|3W|7W)_(left|right|l|r)$", ignore_case = TRUE)) %>%
      str_remove(regex("_(0W|3W|7W)(left|right|l|r)$", ignore_case = TRUE)) %>%
      str_remove(regex("_(EES|EESON|EESOFF|ON|OFF).*$", ignore_case = TRUE))
  )

# 確認
df_long %>%
  count(week)

df_long %>%
  select(raw_id, animal, week, side) %>%
  distinct() %>%
  print(n = 30)

# -----------------------------
# 解析したい12項目
# -----------------------------
metrics_to_test <- c(
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

week_levels <- c("0W", "3W", "7W")

metrics_tbl <- tibble(
  metric_name = metrics_to_test,
  metric_key = str_to_lower(str_squish(metrics_to_test))
)

df_long_12 <- df_long %>%
  mutate(
    metric_key = str_to_lower(str_squish(metric)),
    week = factor(week, levels = week_levels)
  ) %>%
  inner_join(metrics_tbl, by = "metric_key")

# 確認
df_long_12 %>%
  count(metric_name, week)

df_meanLR <- df_long_12 %>%
  group_by(metric_name, animal, week) %>%
  summarise(
    value_mean = ifelse(all(is.na(value)), NA_real_, mean(value, na.rm = TRUE)),
    .groups = "drop"
  )

df_meanLR %>%
  count(metric_name, week)

make_complete_wide <- function(dat) {
  
  wide <- dat %>%
    filter(week %in% week_levels) %>%
    select(animal, week, value_mean) %>%
    mutate(week = as.character(week)) %>%
    pivot_wider(
      names_from = week,
      values_from = value_mean
    )
  
  missing_cols <- setdiff(week_levels, names(wide))
  if (length(missing_cols) > 0) {
    wide[missing_cols] <- NA_real_
  }
  
  wide %>%
    select(animal, all_of(week_levels)) %>%
    drop_na(all_of(week_levels))
}

run_friedman_one <- function(dat) {
  
  wide <- make_complete_wide(dat)
  
  if (nrow(wide) < 2) {
    return(tibble(
      n_animals = nrow(wide),
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_
    ))
  }
  
  res <- friedman.test(
    as.matrix(wide %>% select(all_of(week_levels)))
  )
  
  tibble(
    n_animals = nrow(wide),
    statistic = unname(res$statistic),
    df = unname(res$parameter),
    p_value = res$p.value
  )
}

friedman_results <- df_meanLR %>%
  group_by(metric_name) %>%
  group_modify(~ run_friedman_one(.x)) %>%
  ungroup() %>%
  mutate(
    p_BH_12 = p.adjust(p_value, method = "BH"),
    p_Holm_12 = p.adjust(p_value, method = "holm")
  ) %>%
  arrange(p_value)

friedman_results