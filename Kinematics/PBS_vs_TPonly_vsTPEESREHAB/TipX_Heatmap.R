# =========================================
# Heatmap (ALL animals): L-R tip X over time (first 240 frames)
#  - folder: ~/Desktop/Kinema_raw_data
#  - exclude basename starts with "×"
#  - group = token after first "_": PBS/TP/Comb
#  - Right tip X: col 12, Left tip X: col 32
#  - use first 240 frames (4 sec @ 60 fps)
# Output:
#  - Heatmap_out/Heatmap_LminusR_first240_raw.png
#  - Heatmap_out/Heatmap_LminusR_first240_raw.svg
#  - (optional) Heatmap_out/Heatmap_LminusR_first240_z.png/svg
# =========================================

rm(list = ls())

pkgs <- c("tidyverse", "readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
library(tidyverse)
library(readr)
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")

# -------------------------
# settings
# -------------------------
folder <- "~/Desktop/Kinema_raw_data"
out_dir <- file.path(folder, "Heatmap_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fps <- 60
N_USE <- 240
group_order <- c("PBS", "TP", "Comb")

group_colors <- c(
  "PBS"="#0072B2",
  "TP"="#D55E00",
  "Comb"="#2E8B57"
)

# heatmap mode:
# "raw" = L-R (mean-centered within animal) in mm
# "z"   = z-score within animal (pattern only; amplitude removed)
make_z_version <- TRUE

# -------------------------
# helpers
# -------------------------
extract_group <- function(fname) {
  tok <- str_split(fname, "_", simplify = TRUE)
  if (ncol(tok) < 2) return(NA_character_)
  g <- tok[1, 2]
  case_when(
    str_to_upper(g) == "PBS"  ~ "PBS",
    str_to_upper(g) == "TP"   ~ "TP",
    str_to_lower(g) == "comb" ~ "Comb",
    TRUE ~ NA_character_
  )
}

extract_animal_id <- function(fname) {
  tok <- str_split(fname, "_", simplify = TRUE)
  tok[1, 1]
}

rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))

# -------------------------
# load files
# -------------------------
files <- list.files(folder, pattern="\\.csv$", full.names=TRUE)
files <- files[!str_detect(basename(files), "^×")]
stopifnot(length(files) > 0)

# -------------------------
# per-file: build long data (animal, group, t, value)
# -------------------------
one_file <- function(fp) {
  df <- read_csv(fp, show_col_types = FALSE)
  if (ncol(df) < 32) return(NULL)
  
  df <- df %>% slice_head(n = N_USE)
  
  xR <- suppressWarnings(as.numeric(df[[12]]))
  xL <- suppressWarnings(as.numeric(df[[32]]))
  
  ok <- is.finite(xR) & is.finite(xL)
  if (sum(ok) < N_USE) return(NULL)
  
  sig <- (xL - xR)
  sig_mc <- sig - mean(sig)          # mean-center within animal
  
  t_sec <- (0:(N_USE-1)) / fps
  
  tibble(
    animal = extract_animal_id(basename(fp)),
    group  = extract_group(basename(fp)),
    t = t_sec,
    LminusR_raw = sig_mc
  )
}

dat <- purrr::map_dfr(files, one_file)

# QC
stopifnot(nrow(dat) > 0)
dat <- dat %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, levels = group_order))

# activity for sorting (optional but helpful)
meta <- dat %>%
  group_by(animal, group) %>%
  summarise(
    activity_LminusR = rms(LminusR_raw),
    .groups = "drop"
  )

# -------------------------
# manual order (your requested order)
# -------------------------
manual_order <- c(
  "P8", "P3-1", "P4-1", "P4-2", "P3-2", "B1", "B8", "P6",
  "G3", "G2", "G1", "P1", "B7", "G6", "G5", "P2",
  "B4", "B5", "G7", "G8", "B2"
)

manual_order <- unique(manual_order)  # safety

# order rows: use your manual order, then append any missing animals at the end
present_animals <- meta$animal

missing_in_data <- setdiff(manual_order, present_animals)
if (length(missing_in_data) > 0) {
  message("WARNING: these manual IDs were not found in data: ",
          paste(missing_in_data, collapse = ", "))
}

not_listed <- setdiff(present_animals, manual_order)
if (length(not_listed) > 0) {
  message("NOTE: these animals exist in data but not in manual list; appended at end: ",
          paste(not_listed, collapse = ", "))
}

animal_order <- c(intersect(manual_order, present_animals), not_listed)

dat <- dat %>%
  mutate(animal = factor(animal, levels = animal_order))

# -------------------------
# RAW heatmap (shared FIXED color scale, based on P6 dynamics)
# -------------------------

# FIXED scale (cm) based on your P6 plot
p_lim <- 0.98
lim <- as.numeric(quantile(abs(dat$LminusR_raw), probs = p_lim, na.rm = TRUE))

p_raw <- ggplot(dat, aes(x = t, y = animal, fill = LminusR_raw)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "white", high = "#d7191c",
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish,
    name = "L−R (cm)"
  ) +
  labs(
    title = "Heatmap (all animals): L−R tip X (mean-centered), first 240 frames",
    subtitle = paste0("Fixed color scale: ±", lim, " cm | Sorted by group then activity | fps=", fps,
                      " | 0–", (N_USE-1)/fps, " s"),
    x = "Time (s)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(file.path(out_dir, "Heatmap_LminusR_first240_raw_fixedScale_cm.png"),
       p_raw, width = 9.5, height = 6.0, dpi = 300)
ggsave(file.path(out_dir, "Heatmap_LminusR_first240_raw_fixedScale_cm.svg"),
       p_raw, width = 9.5, height = 6.0)

print(p_raw)