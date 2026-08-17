## ============================================
## STEM121+ fibers at 12, 16, and 20 mm caudal
## TP vs Combination
## Mann–Whitney U test with Holm correction
## ============================================

library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------
## 1) データ入力
## -----------------------------
df <- data.frame(
  distance = c(
    rep("12 mm caudal", 10),
    rep("16 mm caudal", 10),
    rep("20 mm caudal", 10)
  ),
  group = c(
    rep("TP", 5), rep("Combination", 5),
    rep("TP", 5), rep("Combination", 5),
    rep("TP", 5), rep("Combination", 5)
  ),
  animal = c(
    "TP1", "TP2", "TP3", "TP5", "TP6",
    "Comb3", "Comb4", "Comb6", "Comb8", "CombB8",
    "TP1", "TP2", "TP3", "TP5", "TP6",
    "Comb3", "Comb4", "Comb6", "Comb8", "CombB8",
    "TP1", "TP2", "TP3", "TP5", "TP6",
    "Comb3", "Comb4", "Comb6", "Comb8", "CombB8"
  ),
  value = c(
    ## 12 mm caudal
    0.188, 0.080, 0.357, 0.245, 0.053,
    0.441, 0.508, 0.266, 0.389, 0.256,
    
    ## 16 mm caudal
    0.105, 0.177, 0.169, 0.212, 0.113,
    0.199, 0.252, 0.286, 0.254, 0.225,
    
    ## 20 mm caudal
    0.085, 0.086, 0.219, 0.143, 0.123,
    0.236, 0.247, 0.317, 0.246, 0.161
  )
)

df$distance <- factor(
  df$distance,
  levels = c("12 mm caudal", "16 mm caudal", "20 mm caudal")
)

df$group <- factor(
  df$group,
  levels = c("TP", "Combination")
)

## -----------------------------
## 2) 各距離で Mann–Whitney U test
## -----------------------------
test_results <- df %>%
  group_by(distance) %>%
  summarise(
    n_TP = sum(group == "TP"),
    n_Combination = sum(group == "Combination"),
    
    mean_TP = mean(value[group == "TP"]),
    mean_Combination = mean(value[group == "Combination"]),
    
    sem_TP = sd(value[group == "TP"]) / sqrt(n_TP),
    sem_Combination = sd(value[group == "Combination"]) / sqrt(n_Combination),
    
    median_TP = median(value[group == "TP"]),
    median_Combination = median(value[group == "Combination"]),
    
    p_value = wilcox.test(
      value[group == "TP"],
      value[group == "Combination"],
      paired = FALSE,
      exact = TRUE,
      alternative = "two.sided"
    )$p.value,
    
    .groups = "drop"
  )

## -----------------------------
## 3) Holm補正
## -----------------------------
test_results <- test_results %>%
  mutate(
    p_Holm = p.adjust(p_value, method = "holm"),
    significance = case_when(
      p_Holm < 0.001 ~ "***",
      p_Holm < 0.01  ~ "**",
      p_Holm < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

## 結果表示
test_results