# =========================
# Two-way repeated-measures ANOVA for BBB score
# Group × Time
# =========================

install.packages("afex")

# 必要パッケージ
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(afex)

# =========================
# 1. Excel読み込み
# =========================

file_path <- file.choose()

raw <- read_excel(
  file_path,
  sheet = 1,
  col_names = FALSE,
  .name_repair = "unique"
)

names(raw) <- paste0("X", seq_len(ncol(raw)))

# =========================
# 2. Excel blockをlong formatへ変換
# =========================

group_levels <- c("PBS", "TP-only", "Combination")

group_rows <- which(str_trim(as.character(raw$X1)) %in% group_levels)

make_long_block <- function(i) {
  
  start_row <- group_rows[i]
  end_row <- ifelse(i < length(group_rows), group_rows[i + 1] - 1, nrow(raw))
  
  group_name <- str_trim(as.character(raw$X1[start_row]))
  
  # Day header
  day_values <- suppressWarnings(as.numeric(unlist(raw[start_row + 1, -1])))
  valid_cols <- which(!is.na(day_values))
  days <- day_values[valid_cols]
  
  # Data rows
  dat <- raw[(start_row + 2):end_row, c(1, valid_cols + 1)]
  colnames(dat) <- c("Animal", paste0("Day", days))
  
  dat <- dat %>%
    filter(!is.na(Animal), Animal != "") %>%
    mutate(
      Group = group_name,
      Animal_label = as.character(Animal),
      AnimalID = paste(Group, Animal_label, row_number(), sep = "_")
    ) %>%
    pivot_longer(
      cols = starts_with("Day"),
      names_to = "Day",
      values_to = "BBB"
    ) %>%
    mutate(
      Day = as.numeric(str_remove(Day, "Day")),
      BBB = as.numeric(BBB),
      Group = factor(Group, levels = group_levels),
      Day_f = factor(Day)
    )
  
  return(dat)
}

df <- bind_rows(lapply(seq_along(group_rows), make_long_block))

# 確認
table(df$Group)
table(df$Day)
head(df)

# =========================
# 3. Two-way repeated-measures ANOVA
# =========================
# between-subject factor: Group
# within-subject factor: Day
# repeated measure ID: AnimalID

anova_res <- aov_ez(
  id = "AnimalID",
  dv = "BBB",
  data = df,
  between = "Group",
  within = "Day_f",
  type = 3,
  anova_table = list(
    correction = "none",  # 前回の F(28,434) に合わせる
    es = "ges"
  )
)

# 結果表示
anova_res
anova_table <- as.data.frame(anova_res$anova_table)
anova_table$Effect <- rownames(anova_table)

anova_table <- anova_table %>%
  select(Effect, everything())

print(anova_table)

# =========================
# 4. 見やすい形で出力
# =========================

format_p <- function(p) {
  ifelse(p < 0.0001, "<0.0001", signif(p, 3))
}

result_table <- anova_table %>%
  mutate(
    Result = paste0(
      "F(",
      `num Df`, ", ",
      `den Df`, ") = ",
      round(F, 2),
      ", P = ",
      format_p(`Pr(>F)`)
    )
  ) %>%
  select(Effect, Result, `num Df`, `den Df`, F, `Pr(>F)`, ges)

print(result_table)

# =========================
# 5. csv保存
# =========================

write.csv(
  result_table,
  "BBB_two_way_repeated_measures_ANOVA_results.csv",
  row.names = FALSE