library(ggplot2)

# ==========================================================
# 共通関数：p.adjust上位17をGeneRatio順に並べてdotplot
# ==========================================================
make_go_dotplot_gg <- function(ego_obj, title_text, n_terms = 17) {
  
  df <- as.data.frame(ego_obj)
  
  df$GeneRatio_num <- sapply(strsplit(df$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
  # p.adjust上位17を抽出
  df_top <- df[order(df$p.adjust), ][1:n_terms, ]
  
  # その17個をGeneRatio順に並べる
  df_top <- df_top[order(df_top$GeneRatio_num, decreasing = TRUE), ]
  
  # y軸順を固定
  df_top$Description <- factor(
    df_top$Description,
    levels = rev(df_top$Description)
  )
  
  p <- ggplot(
    df_top,
    aes(x = GeneRatio_num, y = Description)
  ) +
    geom_point(
      aes(size = Count, color = p.adjust)
    ) +
    scale_color_gradient(
      low = "#d95f5f",
      high = "#3b82c4",
      name = "p.adjust"
    ) +
    scale_size_continuous(name = "Count") +
    labs(
      title = title_text,
      x = "GeneRatio",
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}

# ==========================================================
# TP
# ==========================================================
p_tp_gg_dot <- make_go_dotplot_gg(
  ego_tp_bp_fc05,
  "GO BP: TP-up genes",
  n_terms = 17
)

print(p_tp_gg_dot)

# ==========================================================
# Combination
# ==========================================================
p_comb_gg_dot <- make_go_dotplot_gg(
  ego_comb_bp_fc05,
  "GO BP: Combination-up genes",
  n_terms = 17
)

print(p_comb_gg_dot)