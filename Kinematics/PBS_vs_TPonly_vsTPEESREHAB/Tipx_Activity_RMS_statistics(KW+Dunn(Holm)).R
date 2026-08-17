# =========================================
# Activity (Tip X) : RMS-based activity from Kinema CSV (FIRST 240 FRAMES)
#  - Use ONLY the first 240 frames (4 sec @ 60 fps)
#  - Exclude files whose basename starts with "×"
#  - Group label = token after first "_" : PBS / TP / Comb
#  - Right tip X: col 12, Left tip X: col 32
# Output:
#  - Activity_out/activity_metrics_all_first240.csv
#  - Activity_out/Activity_boxplot_first240.png
# =========================================

rm(list = ls())

# packages
pkgs <- c("tidyverse", "readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(tidyverse)
library(readr)

# -------------------------
# 0) EDIT HERE
# -------------------------
folder <- "~/Desktop/Kinema_raw_data"
out_dir <- file.path(folder, "Activity_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

group_order <- c("PBS", "TP", "Comb")

group_colors <- c(
  "PBS"  = "#0072B2",
  "TP"   = "#D55E00",
  "Comb" = "#2E8B57"
)

N_USE <- 240  # ★ first 240 frames only

# -------------------------
# 1) Helper functions
# -------------------------
rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))

extract_group <- function(fname) {
  tok <- str_split(fname, "_", simplify = TRUE)
  if (ncol(tok) < 2) return(NA_character_)
  g <- tok[1, 2]
  g <- case_when(
    str_to_upper(g) == "PBS"  ~ "PBS",
    str_to_upper(g) == "TP"   ~ "TP",
    str_to_lower(g) == "comb" ~ "Comb",
    TRUE ~ NA_character_
  )
  g
}

extract_animal_id <- function(fname) {
  tok <- str_split(fname, "_", simplify = TRUE)
  tok[1, 1]
}

# -------------------------
# 2) List CSVs & filter
# -------------------------
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
files <- files[!str_detect(basename(files), "^×")]

stopifnot(length(files) > 0)
message("CSV files found (after excluding ×): ", length(files))

# -------------------------
# 3) Per-file Activity calc (FIRST 240 frames)
# -------------------------
calc_activity_one <- function(file_path) {
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # guard: need enough columns
  if (ncol(df) < 32) {
    return(tibble(
      file = basename(file_path),
      animal = extract_animal_id(basename(file_path)),
      group  = extract_group(basename(file_path)),
      n_frames = nrow(df),
      ok = FALSE,
      activity_rms_L = NA_real_,
      activity_rms_R = NA_real_,
      activity_rms_sumLR = NA_real_,
      reason = "ncol<32"
    ))
  }
  
  # ★ take only first 240 frames
  df240 <- df %>% slice_head(n = N_USE)
  
  xR <- suppressWarnings(as.numeric(df240[[12]]))
  xL <- suppressWarnings(as.numeric(df240[[32]]))
  
  dat <- tibble(xR = xR, xL = xL) %>%
    filter(!is.na(xR) & !is.na(xL))
  
  # ★ require 240 valid frames for fair comparison
  if (nrow(dat) < N_USE) {
    return(tibble(
      file = basename(file_path),
      animal = extract_animal_id(basename(file_path)),
      group  = extract_group(basename(file_path)),
      n_frames = nrow(dat),
      ok = FALSE,
      activity_rms_L = NA_real_,
      activity_rms_R = NA_real_,
      activity_rms_sumLR = NA_real_,
      reason = "less_than_240_valid_frames"
    ))
  }
  
  # mean-center
  xR_mc <- dat$xR - mean(dat$xR, na.rm = TRUE)
  xL_mc <- dat$xL - mean(dat$xL, na.rm = TRUE)
  
  rmsR <- rms(xR_mc)
  rmsL <- rms(xL_mc)
  
  tibble(
    file = basename(file_path),
    animal = extract_animal_id(basename(file_path)),
    group  = extract_group(basename(file_path)),
    n_frames = N_USE,
    ok = TRUE,
    activity_rms_L = rmsL,
    activity_rms_R = rmsR,
    activity_rms_sumLR = rmsL + rmsR,
    reason = NA_character_
  )
}

activity_all <- map_dfr(files, calc_activity_one) %>%
  mutate(group = factor(group, levels = group_order))

# quick check
print(activity_all %>% count(group, ok))
print(activity_all %>% filter(is.na(group)) %>% select(file) %>% head())

# save table
write_csv(activity_all, file.path(out_dir, "activity_metrics_all_first240.csv"))

# -------------------------
# 4) Plot (box + points) : SUM(L+R)
# -------------------------
p_activity <- activity_all %>%
  filter(ok, !is.na(group)) %>%
  ggplot(aes(x = group, y = activity_rms_sumLR, fill = group)) +
  geom_boxplot(
    color = "black",
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.15,
    size = 2.4,
    alpha = 0.9
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Activity (Tip X RMS, L+R) — first 240 frames",
    x = "Group",
    y = "Activity (RMS_L + RMS_R, mean-centered Tip X)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black")
  )

ggsave(
  filename = file.path(out_dir, "Activity_boxplot_first240.png"),
  plot = p_activity,
  width = 5.5, height = 4.2, dpi = 300
)

print(p_activity)

# =========================================
# 5) Statistics: Kruskal-Wallis + Dunn (BH)
# =========================================

# packages for nonparametric stats
if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
library(rstatix)

dat_stat <- activity_all %>%
  filter(ok, !is.na(group)) %>%
  transmute(
    animal,
    group,
    activity = activity_rms_sumLR
  ) %>%
  mutate(group = factor(group, levels = group_order))

stopifnot(nrow(dat_stat) > 0)

# ---- Overall test (3-group): Kruskal-Wallis
kw <- kruskal_test(dat_stat, activity ~ group)

# effect size (epsilon^2; nonparametric analogue)
eff <- kruskal_effsize(dat_stat, activity ~ group)

kw_out <- kw %>%
  left_join(eff, by = character()) %>%   # combine side-by-side
  mutate(
    metric = "Activity_RMS_sumLR",
    test = "Kruskal-Wallis",
    p.value = p,
    p.adj = NA_real_
  ) %>%
  select(metric, test, df, statistic, p.value, effsize, magnitude)

print(kw_out)

write_csv(kw_out, file.path(out_dir, "activity_KW_overall.csv"))

# ---- Post-hoc: Dunn test + BH correction (pairwise)
posthoc <- dat_stat %>%
  dunn_test(activity ~ group, p.adjust.method = "BH") %>%
  arrange(p.adj) %>%
  mutate(metric = "Activity_RMS_sumLR") %>%
  select(metric, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)

print(posthoc)

write_csv(posthoc, file.path(out_dir, "activity_posthoc_Dunn_BH.csv"))

message("Stats saved to: ", out_dir)

# =========================================
# 5) Statistics: Kruskal-Wallis + Dunn (Holm)
# =========================================

# packages for nonparametric stats
if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
library(rstatix)

dat_stat <- activity_all %>%
  filter(ok, !is.na(group)) %>%
  transmute(
    animal,
    group,
    activity = activity_rms_sumLR
  ) %>%
  mutate(group = factor(group, levels = group_order))

stopifnot(nrow(dat_stat) > 0)

# ---- Overall test (3-group): Kruskal-Wallis
kw <- kruskal_test(dat_stat, activity ~ group)

# effect size (epsilon^2; nonparametric analogue)
eff <- kruskal_effsize(dat_stat, activity ~ group)

kw_out <- kw %>%
  left_join(eff, by = character()) %>%   # combine side-by-side
  mutate(
    metric = "Activity_RMS_sumLR",
    test = "Kruskal-Wallis",
    p.value = p,
    p.adj = NA_real_
  ) %>%
  select(metric, test, df, statistic, p.value, effsize, magnitude)

print(kw_out)

write_csv(kw_out, file.path(out_dir, "activity_KW_overall.csv"))

# ---- Post-hoc: Dunn test + Holm correction (pairwise)
posthoc <- dat_stat %>%
  dunn_test(activity ~ group, p.adjust.method = "holm") %>%
  arrange(p.adj) %>%
  mutate(metric = "Activity_RMS_sumLR") %>%
  select(metric, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)

print(posthoc)

write_csv(posthoc, file.path(out_dir, "activity_posthoc_Dunn_BH.csv"))

message("Stats saved to: ", out_dir)