# =========================
# DESeq2::plotPCA for TP-Caud vs Comb-Caud
# =========================

# 必要なら最初にインストール
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("ggplot2")

library(DESeq2)
library(ggplot2)

# 1. 作業フォルダ
# -------------------------
setwd("~/Desktop/bulkRNA_seq")
list.files()

# -------------------------
# 2. データ読み込み
# raw count を含む all.TPM_anno.xls を使う
# -------------------------
dat <- read.delim(
  "all.TPM_anno.xls",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# gene_id を行名に
rownames(dat) <- dat$GeneId

# -------------------------
# 3. TP-Caud / Comb-Caud の raw count 列だけ抽出
# -------------------------
count_cols <- c(
  "TP1-Caud_count",
  "TP2-Caud_count",
  "TP3-Caud_count",
  "TP4-Caud_count",
  "Comb1-Caud_count",
  "Comb2-Caud_count",
  "Comb3-Caud_count",
  "Comb4-Caud_count"
)

count_cols %in% colnames(dat)

count_mat <- dat[, count_cols]

# 数値化
count_mat <- as.data.frame(lapply(count_mat, as.numeric))
rownames(count_mat) <- rownames(dat)

# 整数化（DESeq2用）
count_mat <- round(as.matrix(count_mat))

# -------------------------
# 4. サンプル情報
# -------------------------
coldata <- data.frame(
  row.names = colnames(count_mat),
  group = factor(c("TP", "TP", "TP", "TP",
                   "Combination", "Combination", "Combination", "Combination"))
)

coldata


# -------------------------
# 5. DESeq2 オブジェクト作成
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ group
)

# 低カウント遺伝子を軽く除外
# ここはPCAの安定化のための最低限の前処理
dds <- dds[rowSums(counts(dds)) > 1, ]

# -------------------------
# 6. vst変換
# blind = TRUE はサンプル関係の可視化向け
# -------------------------
vsd <- vst(dds, blind = TRUE)

# -------------------------
# 7. DESeq2標準の plotPCA
# ※ デフォルトでは上位500可変遺伝子
# -------------------------
p1 <- plotPCA(vsd, intgroup = "group") +
  scale_color_manual(values = c(
    "TP" = "#D55E00",
    "Combination" = "#2E8B57"
  )) +
  labs(color = NULL,
       title = "DESeq2 plotPCA: TP-Caud vs Combination-Caud") +
  theme_classic(base_size = 14)

print(p1)

# -------------------------
# 8. PCA座標を取得
# -------------------------
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

print(pcaData)

# PCA座標を保存
write.csv(
  pcaData,
  file = "~/Desktop/DESeq2_plotPCA_scores_TP-Caud_vs_Combination-Caud.csv",
  row.names = TRUE
)

# -------------------------
# 9. 全体の統計（PERMANOVA）
# plotPCA と同じ top 500 variable genes を使う
# -------------------------
# 必要なら最初にインストール
# install.packages("vegan")

library(vegan)

rv <- apply(assay(vsd), 1, var)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]

mat_for_perm <- t(assay(vsd)[select, ])

set.seed(1)

perm_res <- adonis2(
  dist(mat_for_perm) ~ group,
  data = coldata,
  permutations = 9999,
  method = "euclidean"
)

print(perm_res)

capture.output(
  perm_res,
  file = "~/Desktop/DESeq2_plotPCA_PERMANOVA_TP-Caud_vs_Combination-Caud.txt"
)

