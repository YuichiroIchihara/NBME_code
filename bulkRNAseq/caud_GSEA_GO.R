## ============================================
## TP-caud vs Comb-caud
## GSEA (GO Biological Process)
## preranked by DESeq2 stat
##
## このDE表では確認済み：
##   stat > 0 / log2FoldChange > 0  = Comb-caud 高値
## よって
##   NES > 0 = Comb-caud enriched
##   NES < 0 = TP-caud enriched
## ============================================

## -----------------------------
## 0) 必要パッケージ
## -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkgs <- c("clusterProfiler", "enrichplot", "org.Rn.eg.db", "AnnotationDbi", "ggplot2", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)

folder <- "~/Desktop/bulkRNA_seq"

## -----------------------------
## 1) DEファイル読み込み
## -----------------------------
de <- read.delim(
  file.path(folder, "Group_TP-Caud-VS-Comb-Caud_DE_anno.xls"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## ENSEMBL version suffix が付いていても対応できるように除去
## 例: ENSRNOG00000012345.1 -> ENSRNOG00000012345
de$gene_id <- sub("\\..*$", "", rownames(de))

## 念のため必要列確認
print(head(de[, c("gene_id", "GeneSymbol", "stat", "log2FoldChange", "padj")]))

## -----------------------------
## 2) ENSEMBL -> ENTREZID に変換
## -----------------------------
map_df <- AnnotationDbi::select(
  org.Rn.eg.db,
  keys = unique(de$gene_id),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
)

map_df <- map_df[!is.na(map_df$ENTREZID), ]

## DE表と結合
de2 <- merge(
  de,
  map_df,
  by.x = "gene_id",
  by.y = "ENSEMBL",
  all.x = FALSE,
  all.y = FALSE
)

## stat がない/NA のものを除く
de2 <- de2[!is.na(de2$stat), ]

## 同じENTREZIDに複数行ある場合は |stat| 最大のものを採用
de2 <- de2 %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\nMapped genes:", nrow(de2), "\n")

## -----------------------------
## 3) preranked gene list 作成
## -----------------------------
geneList <- de2$stat
names(geneList) <- de2$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

cat("\nSummary of geneList:\n")
print(summary(geneList))

cat("\nTop of geneList:\n")
print(head(geneList))

cat("\nBottom of geneList:\n")
print(tail(geneList))

## -----------------------------
## 4) GSEA: GO Biological Process
## -----------------------------
set.seed(123)

gsea_go <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Rn.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  eps           = 0,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)


## 結果テーブル
gsea_res <- as.data.frame(gsea_go)

## 保存
write.table(
  gsea_res,
  file = file.path(folder, "GSEA_GO_BP_TP-Caud_vs_Comb-Caud.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nTop GSEA results:\n")
print(head(gsea_res[, c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust")], 20))

## -----------------------------
## 5) Comb高値 / TP高値 の上位を確認
## -----------------------------
## この比較では：
## NES > 0  = Comb-caud enriched
## NES < 0  = TP-caud enriched

comb_up <- gsea_res %>%
  filter(NES > 0) %>%
  arrange(p.adjust, desc(NES))

tp_up <- gsea_res %>%
  filter(NES < 0) %>%
  arrange(p.adjust, NES)

cat("\n=== Comb-caud enriched (NES > 0) ===\n")
print(head(comb_up[, c("Description", "NES", "p.adjust")], 20))

cat("\n=== TP-caud enriched (NES < 0) ===\n")
print(head(tp_up[, c("Description", "NES", "p.adjust")], 20))

write.table(
  comb_up,
  file = file.path(folder, "GSEA_GO_BP_Comb_up.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tp_up,
  file = file.path(folder, "GSEA_GO_BP_TP_up.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

make_gsea_dotplot <- function(df, title_text) {
  
  plot_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(NES_abs = abs(NES)) %>%
    arrange(NES_abs)
  
  plot_df$Description <- factor(
    plot_df$Description,
    levels = plot_df$Description
  )
  
  p <- ggplot(plot_df, aes(x = NES_abs, y = Description)) +
    geom_point(
      aes(size = setSize, color = p.adjust)
    ) +
    scale_color_gradientn(
      colors = c("#d95f5f", "#b56576", "#7b6ba8", "#3b82c4"),
      guide = guide_colorbar(reverse = TRUE)
    ) +
    theme_classic(base_size = 11) +
    labs(
      title = title_text,
      x = "|NES|",
      y = NULL,
      size = "Count",
      color = "p.adjust"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  return(p)
}

p_comb_top10 <- make_gsea_dotplot(
  comb_up,
  "Comb-caud enriched GO BP terms"
)

p_tp_top10 <- make_gsea_dotplot(
  tp_up,
  "TP-caud enriched GO BP terms"
)

print(p_comb_top10)
print(p_tp_top10)

## -----------------------------
## 6) Dotplot作成
## term文字なし・Y軸tickあり
## -----------------------------
make_gsea_dotplot <- function(df, title_text) {
  
  plot_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(NES_abs = abs(NES)) %>%
    arrange(NES_abs)
  
  plot_df$Description <- factor(
    plot_df$Description,
    levels = plot_df$Description
  )
  
  p <- ggplot(plot_df, aes(x = NES_abs, y = Description)) +
    geom_point(
      aes(size = setSize, color = p.adjust)
    ) +
    scale_color_gradientn(
      colors = c("#d95f5f", "#b56576", "#7b6ba8", "#3b82c4"),
      guide = guide_colorbar(reverse = TRUE)
    ) +
    theme_classic(base_size = 11) +
    labs(
      title = title_text,
      x = "|NES|",
      y = NULL,
      size = "Count",
      color = "p.adjust"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_blank(),      # term文字だけ消す
      axis.ticks.y = element_line(color = "black"),  # Y軸tickは残す
      axis.text.x = element_text(size = 9, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  return(p)
}

p_comb_top10 <- make_gsea_dotplot(
  comb_up,
  "Comb-caud enriched GO BP terms"
)

p_tp_top10 <- make_gsea_dotplot(
  tp_up,
  "TP-caud enriched GO BP terms"
)

print(p_comb_top10)
print(p_tp_top10)