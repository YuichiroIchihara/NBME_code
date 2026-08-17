# =========================
# DESeq2::plotPCA for TP-Rost vs Comb-Rost
# =========================

library(DESeq2)
library(ggplot2)

# -------------------------
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
# 3. TP-Rost / Comb-Rost の raw count 列だけ抽出
# -------------------------
count_cols <- c(
  "TP1-Rost_count",
  "TP2-Rost_count",
  "TP3-Rost_count",
  "TP4-Rost_count",
  "Comb1-Rost_count",
  "Comb2-Rost_count",
  "Comb3-Rost_count",
  "Comb4-Rost_count"
)

# 列名確認
print(count_cols %in% colnames(dat))

# 足りない列があれば止める
missing_cols <- count_cols[!count_cols %in% colnames(dat)]
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "以下の列が見つかりません:\n",
      paste(missing_cols, collapse = "\n")
    )
  )
}

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
  group = factor(
    c("TP", "TP", "TP", "TP",
      "Combination", "Combination", "Combination", "Combination"),
    levels = c("TP", "Combination")
  )
)

print(coldata)

# -------------------------
# 5. DESeq2 オブジェクト作成
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ group
)

# 低カウント遺伝子を軽く除外
dds <- dds[rowSums(counts(dds)) > 1, ]

# -------------------------
# 6. vst変換
# blind = TRUE はサンプル関係の可視化向け
# -------------------------
vsd <- vst(dds, blind = TRUE)

# -------------------------
# 7. DESeq2標準の plotPCA
# デフォルトでは上位500可変遺伝子
# -------------------------
p1 <- plotPCA(vsd, intgroup = "group") +
  scale_color_manual(values = c(
    "TP" = "#D55E00",
    "Combination" = "#2E8B57"
  )) +
  labs(
    color = NULL,
    title = "DESeq2 plotPCA: TP-Rost vs Combination-Rost"
  ) +
  theme_classic(base_size = 14)

print(p1)

# =========================
# DESeq2::plotPCA for TP-Rost vs Comb-Rost
# =========================

library(DESeq2)
library(ggplot2)

# -------------------------
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
# 3. TP-Rost / Comb-Rost の raw count 列だけ抽出
# -------------------------
count_cols <- c(
  "TP1-Rost_count",
  "TP2-Rost_count",
  "TP3-Rost_count",
  "TP4-Rost_count",
  "Comb1-Rost_count",
  "Comb2-Rost_count",
  "Comb3-Rost_count",
  "Comb4-Rost_count"
)

# 列名確認
print(count_cols %in% colnames(dat))

# 足りない列があれば止める
missing_cols <- count_cols[!count_cols %in% colnames(dat)]
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "以下の列が見つかりません:\n",
      paste(missing_cols, collapse = "\n")
    )
  )
}

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
  group = factor(
    c("TP", "TP", "TP", "TP",
      "Combination", "Combination", "Combination", "Combination"),
    levels = c("TP", "Combination")
  )
)

print(coldata)

# -------------------------
# 5. DESeq2 オブジェクト作成
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ group
)

# 低カウント遺伝子を軽く除外
dds <- dds[rowSums(counts(dds)) > 1, ]

# -------------------------
# 6. vst変換
# blind = TRUE はサンプル関係の可視化向け
# -------------------------
vsd <- vst(dds, blind = TRUE)

# -------------------------
# 7. DESeq2標準の plotPCA
# デフォルトでは上位500可変遺伝子
# -------------------------
p1 <- plotPCA(vsd, intgroup = "group") +
  scale_color_manual(values = c(
    "TP" = "#D55E00",
    "Combination" = "#2E8B57"
  )) +
  labs(
    color = NULL,
    title = "DESeq2 plotPCA: TP-Rost vs Combination-Rost"
  ) +
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
  file = "~/Desktop/DESeq2_plotPCA_scores_TP-Rost_vs_Combination-Rost.csv",
  row.names = TRUE
)

# PCA plot 自体も保存したい場合
ggsave(
  filename = "~/Desktop/DESeq2_plotPCA_TP-Rost_vs_Combination-Rost.pdf",
  plot = p1,
  width = 5.5,
  height = 5
)

ggsave(
  filename = "~/Desktop/DESeq2_plotPCA_TP-Rost_vs_Combination-Rost.png",
  plot = p1,
  width = 5.5,
  height = 5,
  dpi = 300
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
  file = "~/Desktop/DESeq2_plotPCA_PERMANOVA_TP-Rost_vs_Combination-Rost.txt"
)