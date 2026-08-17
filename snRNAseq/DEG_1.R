############################################################
# TP / Combination から human fraction > 0.90 のヒト核を抽出し、
# 軽いQC後に結合してUMAPを作る
#
# 改訂ポイント:
#  1) gene symbol の match() を使わず、Ensembl ID を直接使って
#     human / rat を判定する
#  2) human_fraction > 0.90 でヒト核を抽出する
#  3) rat_fraction も計算して、混入の確認に使えるようにする
#  4) Seuratの標準的な流れでUMAP・クラスタリング・マーカー解析まで行う
#
# 想定:
#  - 10x filtered_feature_bc_matrix を使用
#  - features.tsv.gz の1列目が feature ID (Ensembl ID)
#  - ヒト遺伝子: ENSG...
#  - ラット遺伝子: ENSRNOG... または ENSRNOT...
############################################################

# ==========================================================
# 0. 必要パッケージ
# ==========================================================
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==========================================================
# 1. データフォルダ
# ----------------------------------------------------------
# 10x の filtered_feature_bc_matrix フォルダを指定する
# ==========================================================
tp_dir   <- "~/Desktop/TP_filtered_feature_bc_matrix"
comb_dir <- "~/Desktop/Combination_filtered_feature_bc_matrix"

# ==========================================================
# 2. 10xデータ読み込み
# ----------------------------------------------------------
# Read10X() でカウント行列を読み込む
# ここで得られる行名は gene symbol のことが多いが、
# species 判定は後で features.tsv.gz の Ensembl ID を使って行う
# ==========================================================
tp_counts   <- Read10X(data.dir = tp_dir)
comb_counts <- Read10X(data.dir = comb_dir)

cat("TP counts dim:\n")
print(dim(tp_counts))

cat("Combination counts dim:\n")
print(dim(comb_counts))

# ==========================================================
# 3. features.tsv.gz を読む
# ----------------------------------------------------------
# 通常:
#   V1 = feature ID (Ensembl ID)
#   V2 = gene symbol
#   V3 = feature type
#
# species 判定は gene symbol ではなく Ensembl ID を使う
# ==========================================================
tp_feat <- read.delim(
  gzfile(file.path(tp_dir, "features.tsv.gz")),
  header = FALSE,
  stringsAsFactors = FALSE
)

comb_feat <- read.delim(
  gzfile(file.path(comb_dir, "features.tsv.gz")),
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(tp_feat)   <- c("feature_id", "gene_symbol", "feature_type")
colnames(comb_feat) <- c("feature_id", "gene_symbol", "feature_type")

cat("TP features head:\n")
print(head(tp_feat))

cat("Combination features head:\n")
print(head(comb_feat))

# ==========================================================
# 4. counts 行列の行と features.tsv の行の対応を確認
# ----------------------------------------------------------
# 10x の通常出力では、matrix.mtx の行順と features.tsv.gz の行順は対応する
# そのため、行数が一致していれば、counts 行列の各行に対して
# features.tsv の feature_id を直接割り当てられる
#
# これにより gene symbol の重複やズレを避けられる
# ==========================================================
stopifnot(nrow(tp_counts) == nrow(tp_feat))
stopifnot(nrow(comb_counts) == nrow(comb_feat))

# ----------------------------------------------------------
# species 判定用に、counts 行列の行へ feature_id を対応させる
# 注意:
#   Seurat用の可読性を考えると gene symbol 行名のままでもよいが、
#   species 判定だけは feature_id ベースで行う
# ----------------------------------------------------------
tp_feature_id_in_counts   <- tp_feat$feature_id
comb_feature_id_in_counts <- comb_feat$feature_id

# ==========================================================
# 5. Ensembl ID で human / rat を判定
# ----------------------------------------------------------
# ヒト: ENSG...
# ラット: ENSRNOG... または ENSRNOT...
#
# multi-species reference を想定している
# ==========================================================
tp_is_human <- grepl("^ENSG", tp_feature_id_in_counts)
comb_is_human <- grepl("^ENSG", comb_feature_id_in_counts)

tp_is_rat <- grepl("^ENSRNOG|^ENSRNOT", tp_feature_id_in_counts)
comb_is_rat <- grepl("^ENSRNOG|^ENSRNOT", comb_feature_id_in_counts)

cat("TP human gene rows:\n")
print(table(tp_is_human, useNA = "ifany"))

cat("TP rat gene rows:\n")
print(table(tp_is_rat, useNA = "ifany"))

cat("Combination human gene rows:\n")
print(table(comb_is_human, useNA = "ifany"))

cat("Combination rat gene rows:\n")
print(table(comb_is_rat, useNA = "ifany"))

# ==========================================================
# 6. human fraction / rat fraction を計算
# ----------------------------------------------------------
# human fraction =
#   (ヒト遺伝子にマップされたUMI総数) / (全UMI総数)
#
# rat fraction も同様に計算しておくと、
# 宿主由来RNA混入や mixed nucleus の確認に役立つ
# ==========================================================
tp_total_counts   <- Matrix::colSums(tp_counts)
comb_total_counts <- Matrix::colSums(comb_counts)

tp_human_counts <- Matrix::colSums(
  tp_counts[tp_is_human, , drop = FALSE]
)

comb_human_counts <- Matrix::colSums(
  comb_counts[comb_is_human, , drop = FALSE]
)

tp_rat_counts <- Matrix::colSums(
  tp_counts[tp_is_rat, , drop = FALSE]
)

comb_rat_counts <- Matrix::colSums(
  comb_counts[comb_is_rat, , drop = FALSE]
)

tp_human_fraction   <- tp_human_counts / tp_total_counts
comb_human_fraction <- comb_human_counts / comb_total_counts

tp_rat_fraction   <- tp_rat_counts / tp_total_counts
comb_rat_fraction <- comb_rat_counts / comb_total_counts

cat("TP human_fraction summary:\n")
print(summary(tp_human_fraction))

cat("TP rat_fraction summary:\n")
print(summary(tp_rat_fraction))

cat("Combination human_fraction summary:\n")
print(summary(comb_human_fraction))

cat("Combination rat_fraction summary:\n")
print(summary(comb_rat_fraction))

# ==========================================================
# 7. Seurat object 作成
# ----------------------------------------------------------
# ここではまず緩めの条件で object を作り、
# 後段の subset で human fraction とQC条件をかける
# ==========================================================
tp <- CreateSeuratObject(
  counts = tp_counts,
  project = "TP",
  min.cells = 3,
  min.features = 100
)

comb <- CreateSeuratObject(
  counts = comb_counts,
  project = "Combination",
  min.cells = 3,
  min.features = 100
)

tp$group <- "TP"
comb$group <- "Combination"

# ----------------------------------------------------------
# 先ほど計算した species 比率を metadata に追加
# colnames の順番で対応づける
# ==========================================================
tp$human_fraction <- tp_human_fraction[colnames(tp)]
tp$rat_fraction   <- tp_rat_fraction[colnames(tp)]

comb$human_fraction <- comb_human_fraction[colnames(comb)]
comb$rat_fraction   <- comb_rat_fraction[colnames(comb)]

# ==========================================================
# 8. percent.mt を計算
# ----------------------------------------------------------
# ヒトミトコンドリア遺伝子は一般に "MT-" で始まる
# snRNA-seq では scRNA-seq より percent.mt は低めなことが多いが、
# debris や低品質核の目安として一応見る
# ==========================================================
tp[["percent.mt"]] <- PercentageFeatureSet(tp, pattern = "^MT-")
comb[["percent.mt"]] <- PercentageFeatureSet(comb, pattern = "^MT-")

# ----------------------------------------------------------
# QCの初期確認
# human_fraction / rat_fraction も一緒に見る
# ==========================================================
VlnPlot(
  tp,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "human_fraction", "rat_fraction"),
  ncol = 5,
  pt.size = 0.1
)

VlnPlot(
  comb,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "human_fraction", "rat_fraction"),
  ncol = 5,
  pt.size = 0.1
)

# scatter で species contamination をざっくり確認
FeatureScatter(tp, feature1 = "nCount_RNA", feature2 = "human_fraction")
FeatureScatter(tp, feature1 = "nCount_RNA", feature2 = "rat_fraction")

FeatureScatter(comb, feature1 = "nCount_RNA", feature2 = "human_fraction")
FeatureScatter(comb, feature1 = "nCount_RNA", feature2 = "rat_fraction")

# ==========================================================
# 9. ヒト核抽出 + 軽いQC
# ----------------------------------------------------------
# 今回の条件
#   human_fraction > 0.90
#   nFeature_RNA >= 500
#   percent.mt <= 10
#
# 備考:
#   human_fraction の閾値は 0.90 に緩めているので、
#   後で rat_fraction や marker 発現を見て妥当性確認するとよい
# ==========================================================
tp_human <- subset(
  tp,
  subset = human_fraction > 0.90 &
    nFeature_RNA >= 500 &
    percent.mt <= 10
)

comb_human <- subset(
  comb,
  subset = human_fraction > 0.90 &
    nFeature_RNA >= 500 &
    percent.mt <= 10
)

cat("TP remaining human nuclei:", ncol(tp_human), "\n")
cat("Combination remaining human nuclei:", ncol(comb_human), "\n")

# ----------------------------------------------------------
# 抽出後QC確認
# 抽出後に rat_fraction が高い核が残っていないかも確認するとよい
# ==========================================================
VlnPlot(
  tp_human,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "human_fraction", "rat_fraction"),
  ncol = 5,
  pt.size = 0.1
)

VlnPlot(
  comb_human,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "human_fraction", "rat_fraction"),
  ncol = 5,
  pt.size = 0.1
)

# ==========================================================
# 10. 2群を結合
# ----------------------------------------------------------
# group metadata は merge 後も保持される
# add.cell.ids を付けてセル名の衝突を防ぐ
# ==========================================================
merged_human <- merge(
  x = tp_human,
  y = comb_human,
  add.cell.ids = c("TP", "Combination"),
  project = "HumanNuclei_TP_vs_Combination"
)

cat("Merged object:\n")
print(merged_human)

cat("Cell numbers by group:\n")
print(table(merged_human$group))

# ==========================================================
# 11. 正規化・特徴量選択・スケーリング
# ----------------------------------------------------------
# ここではまず standard workflow を使う
# 比較が主目的なら、後で integration を別途検討してもよい
# ==========================================================
merged_human <- NormalizeData(merged_human)

merged_human <- FindVariableFeatures(
  merged_human,
  selection.method = "vst",
  nfeatures = 2000
)

# 可変遺伝子を中心にスケーリングする方が軽く、一般的
merged_human <- ScaleData(
  merged_human,
  features = VariableFeatures(merged_human)
)

# ==========================================================
# 12. PCA / UMAP / clustering
# ----------------------------------------------------------
# 最初は dims = 1:20 で進める
# 必要に応じて ElbowPlot や PC loadings を見て調整する
# ==========================================================
merged_human <- RunPCA(
  merged_human,
  features = VariableFeatures(merged_human)
)

ElbowPlot(merged_human)

merged_human <- FindNeighbors(merged_human, dims = 1:20)
merged_human <- FindClusters(merged_human, resolution = 0.3)
merged_human <- RunUMAP(merged_human, dims = 1:20)

# ==========================================================
# 13. UMAP描画
# ----------------------------------------------------------
# groupごとの分布と cluster の分かれ方を確認
# ==========================================================
p1 <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.5
) + ggtitle("UMAP by group")

p2 <- DimPlot(
  merged_human,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
) + ggtitle("UMAP by cluster")

p1 + p2

DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.5
)

# ==========================================================
# 14. marker 候補をざっくり確認
# ----------------------------------------------------------
# 移植 iPS 由来神経系細胞で想定される marker を可視化
# cluster annotation の足がかりにする
# ==========================================================
FeaturePlot(merged_human, features = c("DCX", "STMN2", "TUBB3", "MAP2"), ncol = 2)
FeaturePlot(merged_human, features = c("SOX2", "NES", "ELAVL3", "GAP43"), ncol = 2)

FeaturePlot(merged_human, features = c("DCX", "STMN2", "ELAVL3", "MAP2"), ncol = 2)
FeaturePlot(merged_human, features = c("SOX2", "NES", "VIM", "HES1"), ncol = 2)
FeaturePlot(merged_human, features = c("SLC1A3", "GFAP", "AQP4", "OLIG2", "PDGFRA"), ncol = 3)

# 追加で増殖細胞も見ると解釈しやすい
FeaturePlot(merged_human, features = c("MKI67", "TOP2A"), ncol = 2)

# ==========================================================
# 15. RNA assay / layer を確認
# ----------------------------------------------------------
# Seurat v5 では merge 後に layer が分かれていることがある
# marker 解析前に JoinLayers() しておくと安全
# ==========================================================
DefaultAssay(merged_human) <- "RNA"

cat("Layers before join:\n")
print(Layers(merged_human[["RNA"]]))

merged_human[["RNA"]] <- JoinLayers(merged_human[["RNA"]])

cat("Layers after join:\n")
print(Layers(merged_human[["RNA"]]))

# 念のため cluster を Idents に設定
Idents(merged_human) <- "seurat_clusters"

# ==========================================================
# 16. 全clusterのマーカー抽出
# ----------------------------------------------------------
# only.pos = TRUE で各clusterの正マーカーのみ取得
# 閾値はやや緩めにして候補を広めに拾う
# ==========================================================
markers_all <- FindAllMarkers(
  object = merged_human,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10,
  test.use = "wilcox"
)

cat("markers_all dim:\n")
print(dim(markers_all))

cat("markers_all head:\n")
print(head(markers_all))

# ==========================================================
# 17. 各clusterの top marker を抽出
# ----------------------------------------------------------
# avg_log2FC の高い順に top10 を抽出
# heatmap 用にも使う
# ==========================================================
top10_markers <- markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("Top10 markers per cluster:\n")
print(top10_markers)

top10_markers %>%
  select(cluster, gene, avg_log2FC, p_val_adj) %>%
  print(n = 100)

# ==========================================================
# 18. heatmap
# ----------------------------------------------------------
# top marker 遺伝子を使ってクラスタごとの特徴を可視化
# cluster 数が多い場合は top5 に減らしてもよい
# ==========================================================
heatmap_genes <- unique(top10_markers$gene)

cat("Number of heatmap genes:", length(heatmap_genes), "\n")
print(heatmap_genes)

DoHeatmap(
  object = merged_human,
  features = heatmap_genes,
  group.by = "seurat_clusters",
  size = 3
) + NoLegend()



############################################################
# Broad annotation を付与して UMAP を描く
# 対象: merged_human (Seurat object)
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

############################################################
# Broad annotation を metadata に追加
############################################################

# 念のため cluster を文字列化
clusters_chr <- as.character(merged_human$seurat_clusters)

# cluster -> broad class 対応
cluster_to_broad <- c(
  "0"  = "Astroglial/progenitor-like",
  "1"  = "Neuron-like",
  "2"  = "Astroglial/progenitor-like",
  "3"  = "FP/ventral progenitor-like",
  "4"  = "OPC/oligodendroglial-like",
  "5"  = "Neuron-like",
  "6"  = "Non-neural-like",
  "7"  = "Neuron-like",
  "8"  = "Neuron-like",
  "9"  = "Neuron-like",
  "10" = "Neuron-like",
  "11" = "Proliferating cells",
  "12" = "OPC/oligodendroglial-like"
)

# ベクトル化して名前を落とす
broad_vec <- unname(cluster_to_broad[clusters_chr])

# 念のため確認
table(clusters_chr, broad_vec, useNA = "ifany")
sum(is.na(broad_vec))

# metadata に追加
merged_human$broad_class <- broad_vec

# factor順を指定
merged_human$broad_class <- factor(
  merged_human$broad_class,
  levels = c(
    "Neuron-like",
    "Astroglial/progenitor-like",
    "OPC/oligodendroglial-like",
    "FP/ventral progenitor-like",
    "Proliferating cells",
    "Non-neural-like"
  )
)

# 確認
table(merged_human$broad_class, useNA = "ifany")

DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
)

############################################################
# 13クラスターを6つの broad class にまとめ、
# さらに 0-5 の番号で表示するUMAPを作る
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(merged_human) <- "RNA"

# 元クラスター番号を文字列化
merged_human$cluster_num_orig <- as.character(merged_human$seurat_clusters)

# ==========================================================
# 1. 13クラスター -> 6大分類 の対応
# ----------------------------------------------------------
# 必要に応じてここは後で調整してください
# ==========================================================
cluster_to_broad <- c(
  "0"  = "Astroglial/progenitor-like",
  "1"  = "Neuron-like",
  "2"  = "Astroglial/progenitor-like",
  "3"  = "FP/ventral progenitor-like",
  "4"  = "OPC/oligodendroglial-like",
  "5"  = "Neuron-like",
  "6"  = "Non-neural-like",
  "7"  = "Neuron-like",
  "8"  = "Neuron-like",
  "9"  = "Neuron-like",
  "10" = "Neuron-like",
  "11" = "Proliferating cells",
  "12" = "OPC/oligodendroglial-like"
)

# broad class を metadata に追加
merged_human$broad_class <- unname(
  cluster_to_broad[merged_human$cluster_num_orig]
)

# 確認
table(merged_human$cluster_num_orig, merged_human$broad_class, useNA = "ifany")

# ==========================================================
# 2. broad class に 0-5 の番号を振る
# ----------------------------------------------------------
# 好きな順番にできます
# ==========================================================
broad_class_to_num <- c(
  "Neuron-like"                    = "0",
  "Astroglial/progenitor-like"     = "1",
  "OPC/oligodendroglial-like"      = "2",
  "FP/ventral progenitor-like"     = "3",
  "Proliferating cells"            = "4",
  "Non-neural-like"                = "5"
)

merged_human$broad_class_num <- unname(
  broad_class_to_num[merged_human$broad_class]
)

# factor化して順序固定
merged_human$broad_class_num <- factor(
  merged_human$broad_class_num,
  levels = c("0", "1", "2", "3", "4", "5")
)

# 確認
table(merged_human$broad_class, merged_human$broad_class_num, useNA = "ifany")

# ==========================================================
# 3. UMAP by 6-class number
# ==========================================================
Idents(merged_human) <- "broad_class_num"

p_umap_broad_num <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class_num",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.7
) +
  ggtitle("UMAP by broad class number") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad_num)

ggsave(
  "UMAP_by_broad_class_number.pdf",
  plot = p_umap_broad_num,
  width = 7,
  height = 5
)

ggsave(
  "UMAP_by_broad_class_number.png",
  plot = p_umap_broad_num,
  width = 7,
  height = 5,
  dpi = 300
)


# ==========================================================
# 5. 6分類の selected marker genes defining each cluster
#    の DotPlot を作成
# ----------------------------------------------------------
# broad_class_num:
#   0 = Neuron-like
#   1 = Astroglial/progenitor-like
#   2 = OPC/oligodendroglial-like
#   3 = FP/ventral progenitor-like
#   4 = Proliferating cells
#   5 = Non-neural-like
# ==========================================================

# 6分類ごとの代表マーカー候補
selected_markers_6class <- c(
  # 0 Neuron-like
  "DCX", "STMN2", "ELAVL3", "MAP2",
  
  # 1 Astroglial/progenitor-like
  "SOX2", "VIM", "SLC1A3", "AQP4",
  
  # 2 OPC/oligodendroglial-like
  "PDGFRA", "OLIG2", "SOX10", "MOG",
  
  # 3 FP/ventral progenitor-like
  "FOXA1", "FOXA2", "NKX6-1", "ARX",
  
  # 4 Proliferating cells
  "MKI67", "TOP2A", "DLGAP5", "AURKA",
  
  # 5 Non-neural-like
  "FLT1", "KDR", "COL1A2", "CLDN11"
)

# 実際にオブジェクト内に存在する遺伝子だけ残す
selected_markers_6class_present <- selected_markers_6class[
  selected_markers_6class %in% rownames(merged_human)
]

cat("Markers found in merged_human:\n")
print(selected_markers_6class_present)

cat("Markers NOT found in merged_human:\n")
print(setdiff(selected_markers_6class, selected_markers_6class_present))

# broad_class_num を Idents にしておく
Idents(merged_human) <- "broad_class_num"

# DotPlot 作成
p_dot_6class <- DotPlot(
  merged_human,
  features = selected_markers_6class_present,
  group.by = "broad_class_num",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Selected marker genes defining each broad class") +
  xlab("") +
  ylab("Cluster") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_6class)

ggsave(
  "DotPlot_selected_marker_genes_by_broad_class.pdf",
  plot = p_dot_6class,
  width = 10,
  height = 5
)

ggsave(
  "DotPlot_selected_marker_genes_by_broad_class.png",
  plot = p_dot_6class,
  width = 10,
  height = 5,
  dpi = 300
)

############################################################
# 6分類 broad class 用 DotPlot
# 使用マーカー:
# DCX, STMN2, SOX2, SLC1A3, AQP4,
# PDGFRA, OLIG2, MOG,
# FOXA1, NKX6-1,
# MKI67, TOP2A,
# FLT1
#
# 前提:
# - merged_human が存在
# - metadata に broad_class_num が入っている
#   0 = Neuron-like
#   1 = Astroglial/progenitor-like
#   2 = OPC/oligodendroglial-like
#   3 = FP/ventral progenitor-like
#   4 = Proliferating cells
#   5 = Non-neural-like
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(merged_human) <- "RNA"
Idents(merged_human) <- "broad_class_num"

# broad_class_num の順序を固定
merged_human$broad_class_num <- factor(
  merged_human$broad_class_num,
  levels = c("0", "1", "2", "3", "4", "5")
)
Idents(merged_human) <- "broad_class_num"

# ==========================================================
# 1. 使用するマーカー遺伝子
# ==========================================================
selected_markers_6class <- c(
  "DCX", "STMN2",
  "SOX2", "SLC1A3", "AQP4",
  "PDGFRA", "OLIG2",
  "FOXA1", "NKX6-1",
  "MKI67",
  "FLT1"
)

# オブジェクト内に存在する遺伝子のみ使用
selected_markers_present <- selected_markers_6class[
  selected_markers_6class %in% rownames(merged_human)
]

cat("Markers found in merged_human:\n")
print(selected_markers_present)

cat("Markers NOT found in merged_human:\n")
print(setdiff(selected_markers_6class, selected_markers_present))

# ==========================================================
# 2. DotPlot 作成
# ==========================================================
p_dot_6class <- DotPlot(
  object = merged_human,
  features = selected_markers_present,
  group.by = "broad_class_num",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Selected marker genes defining each broad class") +
  xlab("") +
  ylab("Broad class") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_6class)

# ==========================================================
# 3. 保存
# ==========================================================
ggsave(
  filename = "DotPlot_selected_marker_genes_by_broad_class.pdf",
  plot = p_dot_6class,
  width = 10,
  height = 5
)

ggsave(
  filename = "DotPlot_selected_marker_genes_by_broad_class.png",
  plot = p_dot_6class,
  width = 10,
  height = 5,
  dpi = 300
)

# ==========================================================
# 4. broad class 番号と意味の対応表
# ==========================================================
broad_class_table <- data.frame(
  broad_class_num = c("0", "1", "2", "3", "4", "5"),
  broad_class = c(
    "Neuron-like",
    "Astroglial/progenitor-like",
    "OPC/oligodendroglial-like",
    "FP/ventral progenitor-like",
    "Proliferating cells",
    "Non-neural-like"
  )
)

print(broad_class_table)

write.csv(
  broad_class_table,
  file = "broad_class_number_table.csv",
  row.names = FALSE
)




############################################################
# Final broad class annotation + UMAP + DotPlot
#
# broad_class_num:
#   0 = Neuron-like
#   1 = Astroglial/progenitor-like
#   2 = OPC-like
#   3 = FP/ventral progenitor-like
#   4 = Oligodendrocyte-like
#   5 = Non-neural-like
#
# marker genes for DotPlot:
#   Neuron-like                   : DCX, STMN2
#   Astroglial/progenitor-like    : SOX2, SLC1A3, AQP4
#   OPC-like                      : PDGFRA, OLIG2
#   FP/ventral progenitor-like    : FOXA1, NKX6-1
#   Oligodendrocyte-like          : SOX10, MOG, CLDN11
#   Non-neural-like               : FLT1
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(merged_human) <- "RNA"

# 元クラスター番号を文字列化
merged_human$cluster_num_orig <- as.character(merged_human$seurat_clusters)

# ==========================================================
# 1. 13クラスター -> 6 broad class の対応
# ----------------------------------------------------------
# いまの注釈に基づく対応
# ==========================================================
cluster_to_broad <- c(
  "0"  = "Astroglial/progenitor-like",
  "1"  = "Neuron-like",
  "2"  = "Astroglial/progenitor-like",
  "3"  = "FP/ventral progenitor-like",
  "4"  = "OPC-like",
  "5"  = "Neuron-like",
  "6"  = "Non-neural-like",
  "7"  = "Neuron-like",
  "8"  = "Neuron-like",
  "9"  = "Neuron-like",
  "10" = "Neuron-like",
  "11" = "Oligodendrocyte-like",
  "12" = "OPC-like"
)

merged_human$broad_class <- unname(
  cluster_to_broad[merged_human$cluster_num_orig]
)

cat("=== cluster -> broad class ===\n")
print(table(merged_human$cluster_num_orig, merged_human$broad_class, useNA = "ifany"))

# ==========================================================
# 2. broad class に 0-5 の番号を振る
# ==========================================================
broad_class_to_num <- c(
  "Neuron-like"                 = "0",
  "Astroglial/progenitor-like"  = "1",
  "OPC-like"                    = "2",
  "FP/ventral progenitor-like"  = "3",
  "Oligodendrocyte-like"        = "4",
  "Non-neural-like"             = "5"
)

merged_human$broad_class_num <- unname(
  broad_class_to_num[merged_human$broad_class]
)

merged_human$broad_class_num <- factor(
  merged_human$broad_class_num,
  levels = c("0", "1", "2", "3", "4", "5")
)

merged_human$broad_class <- factor(
  merged_human$broad_class,
  levels = c(
    "Neuron-like",
    "Astroglial/progenitor-like",
    "OPC-like",
    "FP/ventral progenitor-like",
    "Oligodendrocyte-like",
    "Non-neural-like"
  )
)

cat("=== broad class num -> broad class ===\n")
print(table(merged_human$broad_class_num, merged_human$broad_class, useNA = "ifany"))

# ==========================================================
# 3. UMAP by broad class number
# ==========================================================
Idents(merged_human) <- "broad_class_num"

p_umap_broad_num <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class_num",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.7
) +
  ggtitle("UMAP by broad class number") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad_num)

ggsave(
  "UMAP_by_broad_class_number.pdf",
  plot = p_umap_broad_num,
  width = 7,
  height = 5
)

ggsave(
  "UMAP_by_broad_class_number.png",
  plot = p_umap_broad_num,
  width = 7,
  height = 5,
  dpi = 300
)

# ==========================================================
# 4. UMAP by broad class name
# ==========================================================
p_umap_broad_name <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.7
) +
  ggtitle("UMAP by broad class") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad_name)

ggsave(
  "UMAP_by_broad_class_name.pdf",
  plot = p_umap_broad_name,
  width = 9,
  height = 5
)

ggsave(
  "UMAP_by_broad_class_name.png",
  plot = p_umap_broad_name,
  width = 9,
  height = 5,
  dpi = 300
)

# ==========================================================
# 5. group別に split した broad class UMAP
# ==========================================================
p_umap_broad_split <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class",
  split.by = "group",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.7
) +
  ggtitle("UMAP by broad class (split by group)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad_split)

ggsave(
  "UMAP_by_broad_class_split_by_group.pdf",
  plot = p_umap_broad_split,
  width = 12,
  height = 5
)

ggsave(
  "UMAP_by_broad_class_split_by_group.png",
  plot = p_umap_broad_split,
  width = 12,
  height = 5,
  dpi = 300
)

# ==========================================================
# 6. broad class 用 DotPlot
# ==========================================================
selected_markers_6class <- c(
  "DCX", "STMN2",
  "SOX2", "SLC1A3", "AQP4",
  "PDGFRA", "OLIG2",
  "FOXA1", "NKX6-1",
  "SOX10", "MOG", "CLDN11",
  "FLT1"
)

selected_markers_present <- selected_markers_6class[
  selected_markers_6class %in% rownames(merged_human)
]

cat("=== markers found ===\n")
print(selected_markers_present)

cat("=== markers NOT found ===\n")
print(setdiff(selected_markers_6class, selected_markers_present))

Idents(merged_human) <- "broad_class_num"

p_dot_6class <- DotPlot(
  object = merged_human,
  features = selected_markers_present,
  group.by = "broad_class_num",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Selected marker genes defining each broad class") +
  xlab("") +
  ylab("Broad class") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_6class)

ggsave(
  "DotPlot_selected_marker_genes_by_broad_class.pdf",
  plot = p_dot_6class,
  width = 10,
  height = 5
)

ggsave(
  "DotPlot_selected_marker_genes_by_broad_class.png",
  plot = p_dot_6class,
  width = 10,
  height = 5,
  dpi = 300
)

# ==========================================================
# 7. broad class 番号対応表
# ==========================================================
broad_class_table <- data.frame(
  broad_class_num = c("0", "1", "2", "3", "4", "5"),
  broad_class = c(
    "Neuron-like",
    "Astroglial/progenitor-like",
    "OPC-like",
    "FP/ventral progenitor-like",
    "Oligodendrocyte-like",
    "Non-neural-like"
  )
)

print(broad_class_table)

write.csv(
  broad_class_table,
  file = "broad_class_number_table.csv",
  row.names = FALSE
)

# ==========================================================
# 8. 各 broad class の細胞数
# ==========================================================
cat("=== cell numbers by broad class ===\n")
print(table(merged_human$broad_class_num))
print(prop.table(table(merged_human$broad_class_num)))

cat("=== cell numbers by group x broad class ===\n")
print(table(merged_human$group, merged_human$broad_class_num))
print(prop.table(table(merged_human$group, merged_human$broad_class_num), margin = 1))

############################################################
# 13クラスターに名前を付けて UMAP を描く
# 対象: merged_human
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(merged_human) <- "RNA"
Idents(merged_human) <- "seurat_clusters"

# 念のため cluster を文字列化して metadata に保存
merged_human$cluster_num <- as.character(merged_human$seurat_clusters)

# ==========================================================
# 1. cluster番号 -> cluster名 の対応表
# ==========================================================
cluster_to_name <- c(
  "0"  = "Astroglial/progenitor-like",
  "1"  = "Neuron-like",
  "2"  = "Ciliated ependymal/radial glia-like",
  "3"  = "FP/ventral progenitor-like",
  "4"  = "OPC-like",
  "5"  = "Peptidergic neuron-like",
  "6"  = "Non-neural-like/endothelial-like",
  "7"  = "Interneuron-like",
  "8"  = "Serotonergic neuron-like",
  "9"  = "Mature neuron-like subtype",
  "10" = "Proliferating cells",
  "11" = "Oligodendrocyte-like",
  "12" = "OPC-like subtype"
)

# metadata に cluster_name を追加
merged_human$cluster_name <- unname(
  cluster_to_name[merged_human$cluster_num]
)

# 確認
cat("=== cluster_num -> cluster_name ===\n")
print(table(merged_human$cluster_num, merged_human$cluster_name, useNA = "ifany"))

# ==========================================================
# 2. 表示順を指定
# ----------------------------------------------------------
# UMAPのlegend順などを揃える
# ==========================================================
cluster_name_levels <- c(
  "Astroglial/progenitor-like",
  "Neuron-like",
  "Ciliated ependymal/radial glia-like",
  "FP/ventral progenitor-like",
  "OPC-like",
  "Peptidergic neuron-like",
  "Non-neural-like/endothelial-like",
  "Interneuron-like",
  "Serotonergic neuron-like",
  "Mature neuron-like subtype",
  "Proliferating cells",
  "Oligodendrocyte-like",
  "OPC-like subtype"
)

merged_human$cluster_name <- factor(
  merged_human$cluster_name,
  levels = cluster_name_levels
)

# ==========================================================
# 3. 名前付きUMAP
# ==========================================================
Idents(merged_human) <- "cluster_name"

p_umap_named <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "cluster_name",
  label = TRUE,
  repel = TRUE,
  label.size = 3.8,
  pt.size = 0.7
) +
  ggtitle("UMAP by annotated cluster") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_named)

ggsave(
  filename = "UMAP_13clusters_named.pdf",
  plot = p_umap_named,
  width = 11,
  height = 6
)

ggsave(
  filename = "UMAP_13clusters_named.png",
  plot = p_umap_named,
  width = 11,
  height = 6,
  dpi = 300
)

############################################################
# 7 broad classes で UMAP と DotPlot を作成
#
# broad_class7_num:
#   0 = Neuron-like
#   1 = Astroglial/progenitor-like
#   2 = OPC-like
#   3 = FP/ventral progenitor-like
#   4 = Oligodendrocyte-like
#   5 = Proliferating cells
#   6 = Atypical neuronal-like
#
# marker genes:
#   Neuron-like                   : DCX, STMN2
#   Astroglial/progenitor-like    : SOX2, SLC1A3, AQP4
#   OPC-like                      : PDGFRA, OLIG2
#   FP/ventral progenitor-like    : FOXA1, NKX6-1
#   Oligodendrocyte-like          : SOX10, MOG, CLDN11
#   Proliferating cells           : MKI67, TOP2A
#   Atypical neuronal-like        : FLT1
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(merged_human) <- "RNA"
merged_human$cluster_num_orig <- as.character(merged_human$seurat_clusters)

# ==========================================================
# 1. 13クラスター -> 7 broad classes の対応
# ----------------------------------------------------------
# 現時点の解釈に基づく対応
# ==========================================================
cluster_to_broad7 <- c(
  "0"  = "Astroglial/progenitor-like",
  "1"  = "Neuron-like",
  "2"  = "Astroglial/progenitor-like",
  "3"  = "FP/ventral progenitor-like",
  "4"  = "OPC-like",
  "5"  = "Neuron-like",
  "6"  = "Atypical neuronal-like",
  "7"  = "Neuron-like",
  "8"  = "Neuron-like",
  "9"  = "Neuron-like",
  "10" = "Proliferating cells",
  "11" = "Oligodendrocyte-like",
  "12" = "Neuron-like"
)

merged_human$broad_class7 <- unname(
  cluster_to_broad7[merged_human$cluster_num_orig]
)

cat("=== cluster -> broad_class7 ===\n")
print(table(merged_human$cluster_num_orig, merged_human$broad_class7, useNA = "ifany"))

# ==========================================================
# 2. broad class に番号を振る
# ==========================================================
broad7_to_num <- c(
  "Neuron-like"                 = "0",
  "Astroglial/progenitor-like"  = "1",
  "OPC-like"                    = "2",
  "FP/ventral progenitor-like"  = "3",
  "Oligodendrocyte-like"        = "4",
  "Proliferating cells"         = "5",
  "Atypical neuronal-like"      = "6"
)

merged_human$broad_class7_num <- unname(
  broad7_to_num[merged_human$broad_class7]
)

merged_human$broad_class7_num <- factor(
  merged_human$broad_class7_num,
  levels = c("0", "1", "2", "3", "4", "5", "6")
)

merged_human$broad_class7 <- factor(
  merged_human$broad_class7,
  levels = c(
    "Neuron-like",
    "Astroglial/progenitor-like",
    "OPC-like",
    "FP/ventral progenitor-like",
    "Oligodendrocyte-like",
    "Proliferating cells",
    "Atypical neuronal-like"
  )
)

cat("=== broad_class7_num -> broad_class7 ===\n")
print(table(merged_human$broad_class7_num, merged_human$broad_class7, useNA = "ifany"))

# ==========================================================
# 3. UMAP by broad class number
# ==========================================================
Idents(merged_human) <- "broad_class7_num"

p_umap_broad7_num <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class7_num",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.7
) +
  ggtitle("UMAP by 7 broad class number") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad7_num)

ggsave(
  "UMAP_by_7broad_class_number.pdf",
  plot = p_umap_broad7_num,
  width = 8,
  height = 5
)

ggsave(
  "UMAP_by_7broad_class_number.png",
  plot = p_umap_broad7_num,
  width = 8,
  height = 5,
  dpi = 300
)

# ==========================================================
# 4. UMAP by broad class name
# ==========================================================
p_umap_broad7_name <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class7",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.7
) +
  ggtitle("UMAP by 7 broad classes") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad7_name)

ggsave(
  "UMAP_by_7broad_class_name.pdf",
  plot = p_umap_broad7_name,
  width = 10,
  height = 5
)

ggsave(
  "UMAP_by_7broad_class_name.png",
  plot = p_umap_broad7_name,
  width = 10,
  height = 5,
  dpi = 300
)

# ==========================================================
# 5. group別 split UMAP
# ==========================================================
p_umap_broad7_split <- DimPlot(
  merged_human,
  reduction = "umap",
  group.by = "broad_class7",
  split.by = "group",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.7
) +
  ggtitle("UMAP by 7 broad classes (split by group)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_broad7_split)

ggsave(
  "UMAP_by_7broad_class_split_by_group.pdf",
  plot = p_umap_broad7_split,
  width = 13,
  height = 5
)

ggsave(
  "UMAP_by_7broad_class_split_by_group.png",
  plot = p_umap_broad7_split,
  width = 13,
  height = 5,
  dpi = 300
)

# ==========================================================
# 6. DotPlot 用マーカー
# ==========================================================
selected_markers_7class <- c(
  "DCX", "STMN2",
  "SOX2", "SLC1A3", "AQP4",
  "PDGFRA", "OLIG2",
  "FOXA1", "NKX6-1",
  "SOX10", "MOG", "CLDN11",
  "MKI67", "TOP2A",
  "FLT1"
)

selected_markers_7class_present <- selected_markers_7class[
  selected_markers_7class %in% rownames(merged_human)
]

cat("=== markers found ===\n")
print(selected_markers_7class_present)

cat("=== markers NOT found ===\n")
print(setdiff(selected_markers_7class, selected_markers_7class_present))

# ==========================================================
# 7. DotPlot
# ==========================================================
Idents(merged_human) <- "broad_class7_num"

p_dot_7class <- DotPlot(
  object = merged_human,
  features = selected_markers_7class_present,
  group.by = "broad_class7_num",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Selected marker genes defining each broad class") +
  xlab("") +
  ylab("Broad class") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_7class)

ggsave(
  "DotPlot_selected_marker_genes_by_7broad_class.pdf",
  plot = p_dot_7class,
  width = 11,
  height = 5
)

ggsave(
  "DotPlot_selected_marker_genes_by_7broad_class.png",
  plot = p_dot_7class,
  width = 11,
  height = 5,
  dpi = 300
)

# ==========================================================
# 9. 各 broad class の細胞数
# ==========================================================
cat("=== cell numbers by 7 broad classes ===\n")
print(table(merged_human$broad_class7_num))
print(prop.table(table(merged_human$broad_class7_num)))

cat("=== cell numbers by group x 7 broad classes ===\n")
print(table(merged_human$group, merged_human$broad_class7_num))
print(prop.table(table(merged_human$group, merged_human$broad_class7_num), margin = 1))

############################################################
# 7 broad classes の構成比 stacked bar plot
#
# 前提:
# - merged_human$group がある
# - merged_human$broad_class7 がある
#
# broad_class7:
#   Neuron-like
#   Astroglial/progenitor-like
#   OPC-like
#   FP/ventral progenitor-like
#   Oligodendrocyte-like
#   Proliferating cells
#   Atypical neuronal-like
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. broad class の順序を固定
# ==========================================================
merged_human$broad_class7 <- factor(
  merged_human$broad_class7,
  levels = c(
    "Neuron-like",
    "Astroglial/progenitor-like",
    "OPC-like",
    "FP/ventral progenitor-like",
    "Oligodendrocyte-like",
    "Proliferating cells",
    "Atypical neuronal-like"
  )
)

merged_human$group <- factor(
  merged_human$group,
  levels = c("TP", "Combination")
)

# ==========================================================
# 1. 細胞数テーブル
# ==========================================================
cell_count_table <- table(merged_human$group, merged_human$broad_class7)
print(cell_count_table)

# データフレーム化
plot_df <- as.data.frame(cell_count_table)
colnames(plot_df) <- c("group", "broad_class7", "cell_count")

print(plot_df)

# ==========================================================
# 2. 割合を計算
# ----------------------------------------------------------
# groupごとに割合を出す
# ==========================================================
plot_df <- plot_df %>%
  group_by(group) %>%
  mutate(
    fraction = cell_count / sum(cell_count),
    percent = fraction * 100
  ) %>%
  ungroup()

print(plot_df)

# ==========================================================
# 3. 積み上げ棒グラフ（割合）
# ==========================================================
p_stack_prop <- ggplot(plot_df, aes(x = group, y = fraction, fill = broad_class7)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Cell composition of 7 broad classes",
    x = "",
    y = "Fraction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_stack_prop)

ggsave(
  "StackedBar_7broad_classes_fraction.pdf",
  plot = p_stack_prop,
  width = 8,
  height = 5
)

ggsave(
  "StackedBar_7broad_classes_fraction.png",
  plot = p_stack_prop,
  width = 8,
  height = 5,
  dpi = 300
)

# ==========================================================
# 4. 積み上げ棒グラフ（細胞数）
# ==========================================================
p_stack_count <- ggplot(plot_df, aes(x = group, y = cell_count, fill = broad_class7)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    title = "Cell number of 7 broad classes",
    x = "",
    y = "Cell number"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_stack_count)

ggsave(
  "StackedBar_7broad_classes_count.pdf",
  plot = p_stack_count,
  width = 8,
  height = 5
)

ggsave(
  "StackedBar_7broad_classes_count.png",
  plot = p_stack_count,
  width = 8,
  height = 5,
  dpi = 300
)



############################################################
# Core neuronal clusters を抽出して再クラスタリング
# 対象 cluster: 1, 5, 7, 8, 9, 12
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

DefaultAssay(merged_human) <- "RNA"
Idents(merged_human) <- "seurat_clusters"

# ==========================================================
# 1. core neuron clusters を抽出
# ==========================================================
core_neuron_clusters <- c("1", "5", "7", "8", "9", "12")

neuron_sub <- subset(
  merged_human,
  idents = core_neuron_clusters
)

cat("Cell numbers in neuron_sub:\n")
print(table(neuron_sub$seurat_clusters))
print(table(neuron_sub$group))

# ==========================================================
# 2. 念のため RNA layer を join
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

if ("counts.1" %in% Layers(neuron_sub[["RNA"]])) {
  neuron_sub[["RNA"]] <- JoinLayers(neuron_sub[["RNA"]])
}

# ==========================================================
# 3. 再解析
# ----------------------------------------------------------
# いったん標準 workflow でやる
# ==========================================================
neuron_sub <- NormalizeData(neuron_sub)
neuron_sub <- FindVariableFeatures(neuron_sub, selection.method = "vst", nfeatures = 2000)
neuron_sub <- ScaleData(neuron_sub, features = VariableFeatures(neuron_sub))
neuron_sub <- RunPCA(neuron_sub, features = VariableFeatures(neuron_sub))

ElbowPlot(neuron_sub)

# まずは 1:15 くらいから開始
neuron_sub <- FindNeighbors(neuron_sub, dims = 1:15)
neuron_sub <- FindClusters(neuron_sub, resolution = 0.3)
neuron_sub <- RunUMAP(neuron_sub, dims = 1:15)

# ==========================================================
# 4. UMAP
# ==========================================================
p1 <- DimPlot(
  neuron_sub,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.6
) + ggtitle("Neuronal subclusters")

p2 <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.6
) + ggtitle("Neuronal subclusters split by group")

p1 + p2

############################################################
# neuronal subcluster markers
############################################################

DefaultAssay(neuron_sub) <- "RNA"
Idents(neuron_sub) <- "seurat_clusters"

markers_neuron <- FindAllMarkers(
  neuron_sub,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10,
  test.use = "wilcox"
)

dim(markers_neuron)
head(markers_neuron)

top10_neuron <- markers_neuron %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

top10_neuron %>%
  select(cluster, gene, avg_log2FC, p_val_adj) %>%
  print(n = 100)

write.csv(markers_neuron, "neuron_sub_all_markers.csv", row.names = FALSE)
write.csv(top10_neuron, "neuron_sub_top10_markers.csv", row.names = FALSE)

############################################################
# neuronal subclusters の top markers heatmap
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(neuron_sub) <- "RNA"
Idents(neuron_sub) <- "seurat_clusters"

# ----------------------------------------------------------
# fold change列名を自動判定
# ----------------------------------------------------------
fc_col <- if ("avg_log2FC" %in% colnames(markers_neuron)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(markers_neuron)) {
  "avg_logFC"
} else {
  stop("fold change列が見つかりません")
}

# ----------------------------------------------------------
# 各cluster top 5 marker
# ----------------------------------------------------------
top5_neuron <- markers_neuron %>%
  group_by(cluster) %>%
  arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

top5_neuron %>%
  select(cluster, gene, all_of(fc_col), p_val_adj) %>%
  print(n = 100)

heatmap_genes_neuron <- unique(top5_neuron$gene)

cat("Number of heatmap genes:\n")
print(length(heatmap_genes_neuron))
print(heatmap_genes_neuron)

# ----------------------------------------------------------
# Heatmap
# ----------------------------------------------------------
p_heat_neuron <- DoHeatmap(
  object = neuron_sub,
  features = heatmap_genes_neuron,
  group.by = "seurat_clusters",
  size = 3
) + 
  NoLegend() +
  ggtitle("Neuronal subcluster top markers")

print(p_heat_neuron)

ggsave(
  "NeuronSubclusters_top5marker_heatmap.pdf",
  plot = p_heat_neuron,
  width = 12,
  height = 8
)

ggsave(
  "NeuronSubclusters_top5marker_heatmap.png",
  plot = p_heat_neuron,
  width = 12,
  height = 8,
  dpi = 300
)

############################################################
# neuronal subcluster annotation 用 FeaturePlot
############################################################

DefaultAssay(neuron_sub) <- "RNA"

# ----------------------------------------------------------
# 1) general neuronal / immature neuronal
# ----------------------------------------------------------
genes_general_neuron <- c("DCX", "STMN2", "ELAVL3", "MAP2", "GAP43")
genes_general_neuron <- genes_general_neuron[genes_general_neuron %in% rownames(neuron_sub)]

p_general <- FeaturePlot(
  neuron_sub,
  features = genes_general_neuron,
  ncol = 3
) & ggtitle("General / immature neuronal markers")

print(p_general)

# ----------------------------------------------------------
# 2) GABAergic / interneuron-like
# ----------------------------------------------------------
genes_gaba <- c("GAD1", "GAD2", "SLC32A1", "NXPH2", "HTR2C")
genes_gaba <- genes_gaba[genes_gaba %in% rownames(neuron_sub)]

p_gaba <- FeaturePlot(
  neuron_sub,
  features = genes_gaba,
  ncol = 3
) & ggtitle("GABAergic / interneuron-like markers")

print(p_gaba)

# ----------------------------------------------------------
# 3) serotonergic-like
# ----------------------------------------------------------
genes_5ht <- c("TPH2", "FEV", "SLC6A4", "SLC18A2", "DDC", "LMX1B")
genes_5ht <- genes_5ht[genes_5ht %in% rownames(neuron_sub)]

p_5ht <- FeaturePlot(
  neuron_sub,
  features = genes_5ht,
  ncol = 3
) & ggtitle("Serotonergic-like markers")

print(p_5ht)

# ----------------------------------------------------------
# 4) peptidergic / secretory neuronal-like
# ----------------------------------------------------------
genes_pept <- c("ADCYAP1", "SCG2", "TAC1", "FAM163A", "SYT10")
genes_pept <- genes_pept[genes_pept %in% rownames(neuron_sub)]

p_pept <- FeaturePlot(
  neuron_sub,
  features = genes_pept,
  ncol = 3
) & ggtitle("Peptidergic / secretory neuronal-like markers")

print(p_pept)

# ----------------------------------------------------------
# 5) ventral / developmental subtype markers
# ----------------------------------------------------------
genes_dev <- c("SIM1", "VSX2", "NKX2-2", "PAX2", "PAX6", "LHX1")
genes_dev <- genes_dev[genes_dev %in% rownames(neuron_sub)]

p_dev <- FeaturePlot(
  neuron_sub,
  features = genes_dev,
  ncol = 3
) & ggtitle("Developmental / ventral subtype markers")

print(p_dev)

# ----------------------------------------------------------
# 保存
# ----------------------------------------------------------
ggsave("NeuronSubclusters_FeaturePlot_general_neuron.pdf", plot = p_general, width = 10, height = 6)
ggsave("NeuronSubclusters_FeaturePlot_gaba.pdf", plot = p_gaba, width = 10, height = 6)
ggsave("NeuronSubclusters_FeaturePlot_serotonergic.pdf", plot = p_5ht, width = 10, height = 8)
ggsave("NeuronSubclusters_FeaturePlot_peptidergic.pdf", plot = p_pept, width = 10, height = 6)
ggsave("NeuronSubclusters_FeaturePlot_developmental.pdf", plot = p_dev, width = 10, height = 8)

############################################################
# annotation 用 DotPlot
############################################################

anno_genes <- c(
  "DCX", "STMN2", "ELAVL3", "MAP2", "GAP43",
  "GAD1", "GAD2", "SLC32A1", "NXPH2", "HTR2C",
  "TPH2", "FEV", "SLC6A4", "SLC18A2", "DDC", "LMX1B",
  "ADCYAP1", "SCG2", "TAC1", "FAM163A", "SYT10",
  "SIM1", "VSX2", "NKX2-2", "PAX2", "PAX6", "LHX1"
)

anno_genes <- anno_genes[anno_genes %in% rownames(neuron_sub)]
print(anno_genes)

p_dot_anno <- DotPlot(
  neuron_sub,
  features = anno_genes,
  group.by = "seurat_clusters",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Neuronal subcluster annotation markers") +
  theme_classic(base_size = 13)

print(p_dot_anno)

ggsave(
  "NeuronSubclusters_annotation_DotPlot.pdf",
  plot = p_dot_anno,
  width = 14,
  height = 6
)

ggsave(
  "NeuronSubclusters_annotation_DotPlot.png",
  plot = p_dot_anno,
  width = 14,
  height = 6,
  dpi = 300
)

############################################################
# Neuronal subclusters に仮ラベルを付けて UMAP を描く
# 対象: neuron_sub
#
# subcluster labels:
#   0 = Immature interneuron-like
#   1 = Peptidergic neuron-like
#   2 = Serotonergic-like
#   3 = Ventral interneuron-like
#   4 = VSX2+ interneuron-like
#   5 = Neuronal subtype 1
#   6 = Inhibitory interneuron-like
#   7 = GABAergic interneuron-like
#   8 = Neuronal subtype 2
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"
Idents(neuron_sub) <- "seurat_clusters"

# 念のため cluster番号を文字列で保存
neuron_sub$subcluster_num <- as.character(neuron_sub$seurat_clusters)

# ==========================================================
# 1. subcluster番号 -> 仮ラベル の対応
# ==========================================================
subcluster_to_name <- c(
  "0" = "Immature interneuron-like",
  "1" = "Peptidergic neuron-like",
  "2" = "Serotonergic-like",
  "3" = "Ventral interneuron-like",
  "4" = "VSX2+ interneuron-like",
  "5" = "Neuronal subtype 1",
  "6" = "Inhibitory interneuron-like",
  "7" = "GABAergic interneuron-like",
  "8" = "Neuronal subtype 2"
)

# metadata に追加
neuron_sub$subcluster_name <- unname(
  subcluster_to_name[neuron_sub$subcluster_num]
)

# 確認
cat("=== subcluster_num -> subcluster_name ===\n")
print(table(neuron_sub$subcluster_num, neuron_sub$subcluster_name, useNA = "ifany"))

# ==========================================================
# 2. 表示順を指定
# ==========================================================
subcluster_name_levels <- c(
  "Immature interneuron-like",
  "Peptidergic neuron-like",
  "Serotonergic-like",
  "Ventral interneuron-like",
  "VSX2+ interneuron-like",
  "Neuronal subtype 1",
  "Inhibitory interneuron-like",
  "GABAergic interneuron-like",
  "Neuronal subtype 2"
)

neuron_sub$subcluster_name <- factor(
  neuron_sub$subcluster_name,
  levels = subcluster_name_levels
)

# ==========================================================
# 3. 名前付きUMAP
# ==========================================================
Idents(neuron_sub) <- "subcluster_name"

p_umap_named <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "subcluster_name",
  label = TRUE,
  repel = TRUE,
  label.size = 3.8,
  pt.size = 0.7
) +
  ggtitle("UMAP of neuronal subclusters") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_named)

ggsave(
  filename = "NeuronSubclusters_named_UMAP.pdf",
  plot = p_umap_named,
  width = 11,
  height = 6
)

ggsave(
  filename = "NeuronSubclusters_named_UMAP.png",
  plot = p_umap_named,
  width = 11,
  height = 6,
  dpi = 300
)

############################################################
# Neuronal subclusters を番号だけで表示する UMAP
# 対象: neuron_sub
############################################################

library(Seurat)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"
Idents(neuron_sub) <- "seurat_clusters"

# 念のため文字列化
neuron_sub$subcluster_num <- as.character(neuron_sub$seurat_clusters)
neuron_sub$subcluster_num <- factor(
  neuron_sub$subcluster_num,
  levels = sort(unique(as.character(neuron_sub$subcluster_num)))
)

# ==========================================================
# 1. 番号だけのUMAP
# ==========================================================
p_umap_num <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "subcluster_num",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.7
) +
  ggtitle("UMAP of neuronal subclusters") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_num)

ggsave(
  filename = "NeuronSubclusters_number_only_UMAP.pdf",
  plot = p_umap_num,
  width = 7,
  height = 6
)

ggsave(
  filename = "NeuronSubclusters_number_only_UMAP.png",
  plot = p_umap_num,
  width = 7,
  height = 6,
  dpi = 300
)

# ==========================================================
# 2. groupでsplitした番号UMAP
# ==========================================================
p_umap_num_split <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "subcluster_num",
  split.by = "group",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.7
) +
  ggtitle("UMAP of neuronal subclusters (split by group)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_umap_num_split)

ggsave(
  filename = "NeuronSubclusters_number_only_UMAP_split_by_group.pdf",
  plot = p_umap_num_split,
  width = 12,
  height = 6
)

ggsave(
  filename = "NeuronSubclusters_number_only_UMAP_split_by_group.png",
  plot = p_umap_num_split,
  width = 12,
  height = 6,
  dpi = 300
)

############################################################
# Neuronal subclusters の split by group シンプルUMAP
# - 色は group のみ
# - ラベルなし
# - できるだけシンプル
############################################################

library(Seurat)
library(ggplot2)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

# group の順序を固定
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

# ==========================================================
# 1. シンプルな split by group UMAP
# ==========================================================
p_umap_group_simple <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.7
) +
  ggtitle("group") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p_umap_group_simple)

ggsave(
  filename = "NeuronSubclusters_UMAP_split_by_group_simple.pdf",
  plot = p_umap_group_simple,
  width = 12,
  height = 6
)

ggsave(
  filename = "NeuronSubclusters_UMAP_split_by_group_simple.png",
  plot = p_umap_group_simple,
  width = 12,
  height = 6,
  dpi = 300
)

p_umap_group_simple <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.7,
  cols = c("Combination" = "#F8766D", "TP" = "#00BFC4")
) +
  ggtitle("group") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p_umap_group_simple)

p_umap_group_simple_nolegend <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.7,
  cols = c("Combination" = "#F8766D", "TP" = "#00BFC4")
) +
  ggtitle("group") +
  NoLegend() +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p_umap_group_simple_nolegend)

p_umap_group_simple <- DimPlot(
  neuron_sub,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  pt.size = 0.7,
  cols = c("Combination" = "#F8766D", "TP" = "#00BFC4")
) +
  ggtitle("group") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_umap_group_simple)

############################################################
# Neuronal subcluster composition stacked bar plot
# 対象: neuron_sub
#
# 前提:
# - neuron_sub$group がある
# - neuron_sub$seurat_clusters が neuronal subcluster を表す
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

# ==========================================================
# 0. 前準備
# ==========================================================
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

neuron_sub$subcluster_num <- as.character(neuron_sub$seurat_clusters)
neuron_sub$subcluster_num <- factor(
  neuron_sub$subcluster_num,
  levels = sort(unique(as.character(neuron_sub$subcluster_num)))
)


# ==========================================================
# 1. 細胞数テーブル
# ==========================================================
cell_count_table <- table(neuron_sub$group, neuron_sub$subcluster_num)
print(cell_count_table)

plot_df <- as.data.frame(cell_count_table)
colnames(plot_df) <- c("group", "subcluster", "cell_count")

print(plot_df)

# ==========================================================
# 2. groupごとの割合を計算
# ==========================================================
plot_df <- plot_df %>%
  group_by(group) %>%
  mutate(
    fraction = cell_count / sum(cell_count),
    percent = fraction * 100
  ) %>%
  ungroup()

print(plot_df)

# ==========================================================
# 3. 色を指定
# ----------------------------------------------------------
# 好きに変更可
# ==========================================================
subcluster_colors <- c(
  "0" = "#F8766D",
  "1" = "#D89000",
  "2" = "#A3A500",
  "3" = "#39B600",
  "4" = "#00BF7D",
  "5" = "#00BFC4",
  "6" = "#00B0F6",
  "7" = "#9590FF",
  "8" = "#E76BF3"
)

# ==========================================================
# 4. 積み上げ棒グラフ（割合）
# ==========================================================
p_stack_prop <- ggplot(plot_df, aes(x = group, y = fraction, fill = subcluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = subcluster_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Neuronal subcluster composition",
    x = "",
    y = "Fraction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_stack_prop)

ggsave(
  "NeuronSubcluster_stackedbar_fraction.pdf",
  plot = p_stack_prop,
  width = 8,
  height = 5
)

ggsave(
  "NeuronSubcluster_stackedbar_fraction.png",
  plot = p_stack_prop,
  width = 8,
  height = 5,
  dpi = 300
)

# ==========================================================
# 5. 積み上げ棒グラフ（細胞数）
# ==========================================================
p_stack_count <- ggplot(plot_df, aes(x = group, y = cell_count, fill = subcluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = subcluster_colors) +
  labs(
    title = "Neuronal subcluster cell number",
    x = "",
    y = "Cell number"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_stack_count)

ggsave(
  "NeuronSubcluster_stackedbar_count.pdf",
  plot = p_stack_count,
  width = 8,
  height = 5
)

ggsave(
  "NeuronSubcluster_stackedbar_count.png",
  plot = p_stack_count,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# Neuronal subcluster composition stacked bar plot
# 凡例を番号ではなく subcluster 名で表示する版
#
# 対象: neuron_sub
# 前提:
# - neuron_sub$group がある
# - neuron_sub$seurat_clusters がある
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

# ==========================================================
# 0. 前準備
# ==========================================================
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

neuron_sub$subcluster_num <- as.character(neuron_sub$seurat_clusters)

# ----------------------------------------------------------
# subcluster 番号 -> 名前
# ----------------------------------------------------------
subcluster_to_name <- c(
  "0" = "Immature interneuron-like",
  "1" = "Peptidergic neuron-like",
  "2" = "Serotonergic-like",
  "3" = "Ventral interneuron-like",
  "4" = "VSX2+ interneuron-like",
  "5" = "Neuronal subtype 1",
  "6" = "Inhibitory interneuron-like",
  "7" = "GABAergic interneuron-like",
  "8" = "Neuronal subtype 2"
)

neuron_sub$subcluster_name <- unname(
  subcluster_to_name[neuron_sub$subcluster_num]
)

# 表示順を固定
subcluster_name_levels <- c(
  "Immature interneuron-like",
  "Peptidergic neuron-like",
  "Serotonergic-like",
  "Ventral interneuron-like",
  "VSX2+ interneuron-like",
  "Neuronal subtype 1",
  "Inhibitory interneuron-like",
  "GABAergic interneuron-like",
  "Neuronal subtype 2"
)

neuron_sub$subcluster_name <- factor(
  neuron_sub$subcluster_name,
  levels = subcluster_name_levels
)

# ==========================================================
# 1. 細胞数テーブル
# ==========================================================
cell_count_table <- table(neuron_sub$group, neuron_sub$subcluster_name)
print(cell_count_table)

plot_df <- as.data.frame(cell_count_table)
colnames(plot_df) <- c("group", "subcluster_name", "cell_count")

print(plot_df)

# ==========================================================
# 2. 割合を計算
# ==========================================================
plot_df <- plot_df %>%
  group_by(group) %>%
  mutate(
    fraction = cell_count / sum(cell_count),
    percent = fraction * 100
  ) %>%
  ungroup()

print(plot_df)

# ==========================================================
# 3. 色を指定
# ----------------------------------------------------------
# 番号UMAPと近い印象の色
# ==========================================================
subcluster_name_colors <- c(
  "Immature interneuron-like"      = "#F8766D",
  "Peptidergic neuron-like"        = "#D89000",
  "Serotonergic-like"              = "#A3A500",
  "Ventral interneuron-like"       = "#39B600",
  "VSX2+ interneuron-like"         = "#00BF7D",
  "Neuronal subtype 1"             = "#00BFC4",
  "Inhibitory interneuron-like"    = "#00B0F6",
  "GABAergic interneuron-like"     = "#9590FF",
  "Neuronal subtype 2"             = "#E76BF3"
)

# ==========================================================
# 4. 積み上げ棒グラフ（割合）
# ==========================================================
p_stack_prop_named <- ggplot(
  plot_df,
  aes(x = group, y = fraction, fill = subcluster_name)
) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = subcluster_name_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Neuronal subcluster composition",
    x = "",
    y = "Fraction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_stack_prop_named)

ggsave(
  "NeuronSubcluster_stackedbar_fraction_namedLegend.pdf",
  plot = p_stack_prop_named,
  width = 10,
  height = 6
)

ggsave(
  "NeuronSubcluster_stackedbar_fraction_namedLegend.png",
  plot = p_stack_prop_named,
  width = 10,
  height = 6,
  dpi = 300
)

# ==========================================================
# 5. 積み上げ棒グラフ（細胞数）
# ==========================================================
p_stack_count_named <- ggplot(
  plot_df,
  aes(x = group, y = cell_count, fill = subcluster_name)
) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = subcluster_name_colors) +
  labs(
    title = "Neuronal subcluster cell number",
    x = "",
    y = "Cell number"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_stack_count_named)

ggsave(
  "NeuronSubcluster_stackedbar_count_namedLegend.pdf",
  plot = p_stack_count_named,
  width = 10,
  height = 6
)

ggsave(
  "NeuronSubcluster_stackedbar_count_namedLegend.png",
  plot = p_stack_count_named,
  width = 10,
  height = 6,
  dpi = 300
)

############################################################
# Simple DotPlot for neuronal subtypes
############################################################

library(Seurat)
library(ggplot2)

DefaultAssay(neuron_sub) <- "RNA"
Idents(neuron_sub) <- "seurat_clusters"

simple_genes <- c(
  "LHX1",
  "ADCYAP1", "SYT10",
  "TPH2", "SLC6A4",
  "SIM1", "NKX2-2",
  "VSX2",
  "KMO",
  "SLC6A5",
  "HTR2C", "GAD1",
  "NFATC1"
)

simple_genes <- simple_genes[simple_genes %in% rownames(neuron_sub)]
print(simple_genes)

p_dot_simple <- DotPlot(
  neuron_sub,
  features = simple_genes,
  group.by = "seurat_clusters",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Selected marker genes for neuronal subtypes") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1)
  )

print(p_dot_simple)

ggsave(
  "NeuronSubclusters_simple_DotPlot.pdf",
  plot = p_dot_simple,
  width = 9,
  height = 5
)

ggsave(
  "NeuronSubclusters_simple_DotPlot.png",
  plot = p_dot_simple,
  width = 9,
  height = 5,
  dpi = 300
)

############################################################
# TP vs Combination DEG (overall in neuron_sub)
# volcano plot
#
# 比較:
#   ident.1 = Combination
#   ident.2 = TP
#
# 解釈:
#   avg_log2FC > 0  : Combinationで高い
#   avg_log2FC < 0  : TPで高い
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

if (length(Layers(neuron_sub[["RNA"]])) > 1) {
  neuron_sub[["RNA"]] <- JoinLayers(neuron_sub[["RNA"]])
}

Idents(neuron_sub) <- "group"
table(Idents(neuron_sub))

# ==========================================================
# 1. DEG解析
# ==========================================================
deg_comb_vs_tp <- FindMarkers(
  object = neuron_sub,
  assay = "RNA",
  ident.1 = "Combination",
  ident.2 = "TP",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "wilcox"
)

# gene名を列として追加
deg_comb_vs_tp <- deg_comb_vs_tp %>%
  tibble::rownames_to_column(var = "gene")

# fold change列名を自動判定
fc_col <- if ("avg_log2FC" %in% colnames(deg_comb_vs_tp)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(deg_comb_vs_tp)) {
  "avg_logFC"
} else {
  stop("fold change列が見つかりません")
}

write.csv(
  deg_comb_vs_tp,
  file = "DEG_neuronSub_Combination_vs_TP.csv",
  row.names = FALSE
)

head(deg_comb_vs_tp)

cat("DEG result dim:\n")
print(dim(deg_comb_vs_tp))

cat("Top genes:\n")
print(head(deg_comb_vs_tp[order(deg_comb_vs_tp$p_val_adj), ], 20))

# ==========================================================
# 2. volcano plot 用データ作成
# ----------------------------------------------------------
# 閾値:
#   padj < 0.05
#   |log2FC| > 1   (= fold change > 2)
# ==========================================================
padj_cutoff <- 0.05
log2fc_cutoff <- 1

volcano_df <- deg_comb_vs_tp %>%
  mutate(
    p_val_adj_plot = ifelse(is.na(p_val_adj), 1, p_val_adj),
    p_val_adj_plot = ifelse(p_val_adj_plot == 0, 1e-300, p_val_adj_plot),
    neglog10_padj = -log10(p_val_adj_plot),
    significance = case_when(
      p_val_adj_plot < padj_cutoff & .data[[fc_col]] >  log2fc_cutoff ~ "Up in Combination",
      p_val_adj_plot < padj_cutoff & .data[[fc_col]] < -log2fc_cutoff ~ "Up in TP",
      TRUE ~ "Not significant"
    )
  )

table(volcano_df$significance)

# ==========================================================
# 3. ラベルする遺伝子を選ぶ
# ----------------------------------------------------------
# ここでは
#   - Combination側 上位
#   - TP側 上位
# をそれぞれ padj順で数個ずつラベル
#
# 必要なら manually 指定も可能
# ==========================================================
n_label_each_side <- 8

label_genes_comb <- volcano_df %>%
  filter(significance == "Up in Combination") %>%
  arrange(p_val_adj_plot) %>%
  slice_head(n = n_label_each_side)

label_genes_tp <- volcano_df %>%
  filter(significance == "Up in TP") %>%
  arrange(p_val_adj_plot) %>%
  slice_head(n = n_label_each_side)

label_df <- bind_rows(label_genes_comb, label_genes_tp)

print(label_df %>% select(gene, all_of(fc_col), p_val_adj_plot, significance))

# ==========================================================
# 4. volcano plot
# ----------------------------------------------------------
# 色:
#   Combination高発現 = red
#   TP高発現          = blue
#   その他            = grey
# ==========================================================
p_volcano <- ggplot(
  volcano_df,
  aes(x = .data[[fc_col]], y = neglog10_padj)
) +
  geom_point(
    aes(color = significance),
    size = 2.2,
    alpha = 0.8
  ) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
             linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", linewidth = 0.7) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 5,
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "black",
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Up in Combination" = "#D55E00",
      "Up in TP" = "#0072B2",
      "Not significant" = "grey75"
    )
  ) +
  labs(
    title = "Combination vs TP in neuronal cells",
    x = "avg_log2FC",
    y = "-log10 adjusted P value"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p_volcano)

ggsave(
  filename = "Volcano_neuronSub_Combination_vs_TP.pdf",
  plot = p_volcano,
  width = 9,
  height = 8
)

ggsave(
  filename = "Volcano_neuronSub_Combination_vs_TP.png",
  plot = p_volcano,
  width = 9,
  height = 8,
  dpi = 300
)

library(ggplot2)
library(ggrepel)
library(dplyr)

# ==========================================================
# 閾値
# ==========================================================
padj_cutoff <- 0.05
log2fc_cutoff <- 1   # = fold change 2倍

# ==========================================================
# volcano plot 用データ
# ==========================================================
volcano_df <- deg_comb_vs_tp %>%
  mutate(
    p_val_adj_plot = ifelse(is.na(p_val_adj), 1, p_val_adj),
    p_val_adj_plot = ifelse(p_val_adj_plot == 0, 1e-300, p_val_adj_plot),
    neglog10_padj = -log10(p_val_adj_plot),
    significance = case_when(
      p_val_adj_plot < padj_cutoff & .data[[fc_col]] < -log2fc_cutoff ~ "Up in TP",
      p_val_adj_plot < padj_cutoff & .data[[fc_col]] >  log2fc_cutoff ~ "Up in Combination",
      TRUE ~ "Not significant"
    )
  )

# ==========================================================
# TP up 遺伝子を全部ラベル
# ==========================================================
label_df_tp <- volcano_df %>%
  filter(significance == "Up in TP")

# 必要なら Combination側も全部ラベルする場合
# label_df_comb <- volcano_df %>%
#   filter(significance == "Up in Combination")

# 今回は TP up だけ全部ラベル
label_df <- label_df_tp

cat("Number of TP-up labeled genes:\n")
print(nrow(label_df))

# ==========================================================
# volcano plot
# ==========================================================
p_volcano_tpblue <- ggplot(
  volcano_df,
  aes(x = .data[[fc_col]], y = neglog10_padj)
) +
  geom_point(
    aes(color = significance),
    size = 2.4,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 5.5,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "black",
    segment.size = 0.5,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 2
  ) +
  scale_color_manual(
    values = c(
      "Up in TP" = "#1F4FFF",
      "Up in Combination" = "grey75",
      "Not significant" = "grey75"
    )
  ) +
  labs(
    title = "Combination vs TP in neuronal cells",
    x = "avg_log2FC",
    y = "-log10 adjusted P value"
  ) +
  theme_classic(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 26),
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(p_volcano_tpblue)

ggsave(
  "Volcano_neuronSub_Combination_vs_TP_TPblue_alllabels.pdf",
  plot = p_volcano_tpblue,
  width = 12,
  height = 8
)

ggsave(
  "Volcano_neuronSub_Combination_vs_TP_TPblue_alllabels.png",
  plot = p_volcano_tpblue,
  width = 12,
  height = 8,
  dpi = 300
)

############################################################
# raw p < 0.05, log2FC > 0.5 / < -0.5 で gene を抽出し、
# GO Biological Process解析を行う
#
# 前提:
# - deg_comb_vs_tp が存在
# - fc_col が定義済み
#   (avg_log2FC または avg_logFC)
#
# 解釈:
#   avg_log2FC > 0 : Combination で高い
#   avg_log2FC < 0 : TP で高い
############################################################

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ==========================================================
# 0. 前提確認
# ==========================================================
fc_col <- if ("avg_log2FC" %in% colnames(deg_comb_vs_tp)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(deg_comb_vs_tp)) {
  "avg_logFC"
} else {
  stop("fold change列が見つかりません")
}

if (!"p_val" %in% colnames(deg_comb_vs_tp)) {
  stop("deg_comb_vs_tp に p_val 列がありません")
}

# ==========================================================
# 1. gene list 作成
# ==========================================================
genes_comb_fc05 <- deg_comb_vs_tp %>%
  filter(
    p_val < 0.05,
    .data[[fc_col]] > 0.5
  ) %>%
  arrange(desc(.data[[fc_col]])) %>%
  pull(gene)

genes_tp_fc05 <- deg_comb_vs_tp %>%
  filter(
    p_val < 0.05,
    .data[[fc_col]] < -0.5
  ) %>%
  arrange(.data[[fc_col]]) %>%
  pull(gene)

cat("Number of Combination-up genes:\n")
print(length(genes_comb_fc05))

cat("Number of TP-up genes:\n")
print(length(genes_tp_fc05))

# ==========================================================
# 2. human symbol に寄せる
# ==========================================================
genes_comb_fc05_human <- toupper(genes_comb_fc05)
genes_tp_fc05_human   <- toupper(genes_tp_fc05)

# ==========================================================
# 3. SYMBOL -> ENTREZID
# ==========================================================
gene_df_comb_fc05 <- bitr(
  genes_comb_fc05_human,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_df_tp_fc05 <- bitr(
  genes_tp_fc05_human,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

cat("Mapped Combination-up genes:\n")
print(nrow(gene_df_comb_fc05))

cat("Mapped TP-up genes:\n")
print(nrow(gene_df_tp_fc05))

comb_fc05_entrez <- unique(gene_df_comb_fc05$ENTREZID)
tp_fc05_entrez   <- unique(gene_df_tp_fc05$ENTREZID)

# ==========================================================
# 4. enrichGO
# ==========================================================
run_go <- function(gene_ids, ontology = "BP") {
  enrichGO(
    gene          = gene_ids,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.3,
    readable      = TRUE
  )
}

ego_comb_bp_fc05 <- run_go(comb_fc05_entrez, "BP")
ego_tp_bp_fc05   <- run_go(tp_fc05_entrez, "BP")

ego_comb_bp_fc05_df <- as.data.frame(ego_comb_bp_fc05)
ego_tp_bp_fc05_df   <- as.data.frame(ego_tp_bp_fc05)

cat("=== Combination BP ===\n")
print(head(ego_comb_bp_fc05_df, 20))

cat("=== TP BP ===\n")
print(head(ego_tp_bp_fc05_df, 20))

write.csv(ego_comb_bp_fc05_df, "GO_BP_Combination_rawp0.05_log2FC0.5.csv", row.names = FALSE)
write.csv(ego_tp_bp_fc05_df,   "GO_BP_TP_rawp0.05_log2FCm0.5.csv", row.names = FALSE)

# ==========================================================
# 5. dotplot
# ==========================================================
make_go_dotplot <- function(ego_obj, title_text, out_prefix, show_n = 20) {
  df <- as.data.frame(ego_obj)
  
  if (nrow(df) == 0) {
    message("No enriched terms for: ", title_text)
    return(NULL)
  }
  
  p <- dotplot(
    ego_obj,
    showCategory = min(show_n, nrow(df))
  ) +
    ggtitle(title_text) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  print(p)
  
  ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 6, dpi = 300)
  
  return(p)
}

p_comb_bp_fc05 <- make_go_dotplot(
  ego_comb_bp_fc05,
  "GO BP: Combination-up genes (raw p<0.05, log2FC>0.5)",
  "GO_BP_dotplot_Combination_rawp0.05_log2FC0.5"
)

p_tp_bp_fc05 <- make_go_dotplot(
  ego_tp_bp_fc05,
  "GO BP: TP-up genes (raw p<0.05, log2FC<-0.5)",
  "GO_BP_dotplot_TP_rawp0.05_log2FCm0.5"
)

############################################################
# Axon guidance / neurite pathfinding genes
# DotPlot + ModuleScore + boxplot + statistics
#
# 対象:
#   neuron_sub
#
# 色:
#   Combination = #2E8B57
#   TP          = #D55E00
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

# group順を固定
# DotPlotで TP を上、Combination を下にしたいので
# factor順は Combination -> TP にしておく
# （DotPlotのy軸は下から上に並ぶため）
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

# ==========================================================
# 1. geneset定義
# ----------------------------------------------------------
# ユーザー指定:
# DCC, NRP1, NRP2, PLXNA1, ROBO1, ROBO2, EPHA4, NRCAM,
# FAT3, IGF1R, PTN, RTN4RL1
############################################################
axon_guidance_genes <- c(
  "DCC", "NRP1", "NRP2", "PLXNA1",
  "ROBO1", "ROBO2", "EPHA4", "NRCAM",
  "FAT3", "IGF1R", "PTN", "RTN4RL1"
)

# 実際に存在する遺伝子のみ使う
axon_guidance_genes_present <- axon_guidance_genes[
  axon_guidance_genes %in% rownames(neuron_sub)
]

cat("Genes found in neuron_sub:\n")
print(axon_guidance_genes_present)

cat("Genes NOT found in neuron_sub:\n")
print(setdiff(axon_guidance_genes, axon_guidance_genes_present))

# ==========================================================
# 2. DotPlot
# ----------------------------------------------------------
# group別に平均発現と発現割合を可視化
# TP を上、Combination を下に表示
# ==========================================================
p_dot_axon_guidance <- DotPlot(
  object = neuron_sub,
  features = axon_guidance_genes_present,
  group.by = "group",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Axon guidance / neurite pathfinding genes") +
  xlab("") +
  ylab("Group") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_axon_guidance)

ggsave(
  filename = "DotPlot_axon_guidance_pathfinding_genes_by_group.pdf",
  plot = p_dot_axon_guidance,
  width = 10,
  height = 4.5
)

ggsave(
  filename = "DotPlot_axon_guidance_pathfinding_genes_by_group.png",
  plot = p_dot_axon_guidance,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ==========================================================
# 3. ModuleScore を計算
# ----------------------------------------------------------
# AddModuleScore は name="AxonGuidanceScore" の場合、
# metadata列名は通常 "AxonGuidanceScore1" になる
# ==========================================================
neuron_sub <- AddModuleScore(
  object = neuron_sub,
  features = list(axon_guidance_genes_present),
  name = "AxonGuidanceScore"
)

# 列名確認
cat("Metadata columns containing AxonGuidanceScore:\n")
print(grep("AxonGuidanceScore", colnames(neuron_sub@meta.data), value = TRUE))

score_col <- grep("AxonGuidanceScore", colnames(neuron_sub@meta.data), value = TRUE)[1]


# ==========================================================
# 4. scoreの分布確認
# ==========================================================
VlnPlot(
  object = neuron_sub,
  features = score_col,
  group.by = "group",
  pt.size = 0.1
)

# ==========================================================
# 5. boxplot用データフレーム
# ==========================================================
plot_df_score <- neuron_sub@meta.data %>%
  dplyr::select(group, all_of(score_col)) %>%
  dplyr::rename(module_score = all_of(score_col))

head(plot_df_score)

# ==========================================================
# 6. 統計
# ----------------------------------------------------------
# 基本:
#   2群比較なので Wilcoxon rank-sum test
#
# 必要なら t-test も併記
# ==========================================================
wilcox_res <- wilcox.test(
  module_score ~ group,
  data = plot_df_score
)

ttest_res <- t.test(
  module_score ~ group,
  data = plot_df_score
)

cat("=== Wilcoxon test ===\n")
print(wilcox_res)

cat("=== t-test ===\n")
print(ttest_res)

# p値をデータフレーム化して保存
stats_df <- data.frame(
  test = c("Wilcoxon rank-sum", "Welch t-test"),
  p_value = c(wilcox_res$p.value, ttest_res$p.value)
)

print(stats_df)

write.csv(
  stats_df,
  file = "AxonGuidance_ModuleScore_TP_vs_Combination_statistics.csv",
  row.names = FALSE
)

# ==========================================================
# 7. boxplot
# ----------------------------------------------------------
# 色:
#   Combination #2E8B57
#   TP          #D55E00
#
# 図中に統計値は表示しない
# ==========================================================
# Boxplotでは見慣れた順序として TP -> Combination に戻す
plot_df_score$group <- factor(
  plot_df_score$group,
  levels = c("TP", "Combination")
)

group_colors <- c(
  "Combination" = "#2E8B57",
  "TP" = "#D55E00"
)

p_box_score <- ggplot(
  plot_df_score,
  aes(x = group, y = module_score, fill = group)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.18,
    size = 1.2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Axon guidance / neurite pathfinding module score",
    x = "",
    y = "Module score"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_box_score)

ggsave(
  filename = "Boxplot_AxonGuidance_ModuleScore_TP_vs_Combination.pdf",
  plot = p_box_score,
  width = 6,
  height = 5
)

ggsave(
  filename = "Boxplot_AxonGuidance_ModuleScore_TP_vs_Combination.png",
  plot = p_box_score,
  width = 6,
  height = 5,
  dpi = 300
)

# ==========================================================
# 8. 平均・中央値も保存
# ==========================================================
summary_df <- plot_df_score %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(module_score, na.rm = TRUE),
    sd = sd(module_score, na.rm = TRUE),
    median = median(module_score, na.rm = TRUE),
    IQR = IQR(module_score, na.rm = TRUE)
  )

print(summary_df)

write.csv(
  summary_df,
  file = "AxonGuidance_ModuleScore_TP_vs_Combination_summary.csv",
  row.names = FALSE
)


############################################################
# Axon outgrowth / regeneration-associated genes
# DotPlot + ModuleScore + boxplot + statistics
#
# 対象:
#   neuron_sub
#
# 色:
#   Combination = #2E8B57
#   TP          = #D55E00
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

# DotPlotで TP を上、Combination を下にしたいので
# factor順は Combination -> TP
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

# ==========================================================
# 1. geneset定義
# ==========================================================
axon_outgrowth_genes <- c(
  "GAP43", "STMN2", "MAP1B", "DCX", "TUBA1A", "TUBB3"
)

# 実際に存在する遺伝子のみ使う
axon_outgrowth_genes_present <- axon_outgrowth_genes[
  axon_outgrowth_genes %in% rownames(neuron_sub)
]

cat("Genes found in neuron_sub:\n")
print(axon_outgrowth_genes_present)

cat("Genes NOT found in neuron_sub:\n")
print(setdiff(axon_outgrowth_genes, axon_outgrowth_genes_present))

# ==========================================================
# 2. DotPlot
# ----------------------------------------------------------
# 右側凡例:
#   上 = Percent Expression
#   下 = Average Expression
# ==========================================================
p_dot_axon_outgrowth <- DotPlot(
  object = neuron_sub,
  features = axon_outgrowth_genes_present,
  group.by = "group",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Axon outgrowth / regeneration-associated genes") +
  xlab("") +
  ylab("Group") +
  guides(
    size  = guide_legend(order = 1),     # Percent Expression を上
    color = guide_colorbar(order = 2)    # Average Expression を下
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_axon_outgrowth)

ggsave(
  filename = "DotPlot_axon_outgrowth_regeneration_genes_by_group.pdf",
  plot = p_dot_axon_outgrowth,
  width = 8,
  height = 4.5
)

ggsave(
  filename = "DotPlot_axon_outgrowth_regeneration_genes_by_group.png",
  plot = p_dot_axon_outgrowth,
  width = 8,
  height = 4.5,
  dpi = 300
)

# ==========================================================
# 3. ModuleScore を計算
# ==========================================================
neuron_sub <- AddModuleScore(
  object = neuron_sub,
  features = list(axon_outgrowth_genes_present),
  name = "AxonOutgrowthScore"
)

# 列名確認
cat("Metadata columns containing AxonOutgrowthScore:\n")
print(grep("AxonOutgrowthScore", colnames(neuron_sub@meta.data), value = TRUE))

score_col <- grep("AxonOutgrowthScore", colnames(neuron_sub@meta.data), value = TRUE)[1]

# ==========================================================
# 4. scoreの分布確認
# ==========================================================
VlnPlot(
  object = neuron_sub,
  features = score_col,
  group.by = "group",
  pt.size = 0.1
)

# ==========================================================
# 5. boxplot用データフレーム
# ==========================================================
plot_df_score <- neuron_sub@meta.data %>%
  dplyr::select(group, all_of(score_col)) %>%
  dplyr::rename(module_score = all_of(score_col))

head(plot_df_score)

# ==========================================================
# 6. 統計
# ==========================================================
wilcox_res <- wilcox.test(
  module_score ~ group,
  data = plot_df_score
)

ttest_res <- t.test(
  module_score ~ group,
  data = plot_df_score
)

cat("=== Wilcoxon test ===\n")
print(wilcox_res)

cat("=== t-test ===\n")
print(ttest_res)

stats_df <- data.frame(
  test = c("Wilcoxon rank-sum", "Welch t-test"),
  p_value = c(wilcox_res$p.value, ttest_res$p.value)
)

print(stats_df)

write.csv(
  stats_df,
  file = "AxonOutgrowth_ModuleScore_TP_vs_Combination_statistics.csv",
  row.names = FALSE
)

# ==========================================================
# 7. boxplot
# ----------------------------------------------------------
# 図中に統計値は表示しない
# ==========================================================
plot_df_score$group <- factor(
  plot_df_score$group,
  levels = c("TP", "Combination")
)

group_colors <- c(
  "Combination" = "#2E8B57",
  "TP" = "#D55E00"
)

p_box_score <- ggplot(
  plot_df_score,
  aes(x = group, y = module_score, fill = group)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.18,
    size = 1.2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Axon outgrowth / regeneration-associated module score",
    x = "",
    y = "Module score"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_box_score)

ggsave(
  filename = "Boxplot_AxonOutgrowth_ModuleScore_TP_vs_Combination.pdf",
  plot = p_box_score,
  width = 6,
  height = 5
)

ggsave(
  filename = "Boxplot_AxonOutgrowth_ModuleScore_TP_vs_Combination.png",
  plot = p_box_score,
  width = 6,
  height = 5,
  dpi = 300
)

# ==========================================================
# 8. 平均・中央値も保存
# ==========================================================
summary_df <- plot_df_score %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(module_score, na.rm = TRUE),
    sd = sd(module_score, na.rm = TRUE),
    median = median(module_score, na.rm = TRUE),
    IQR = IQR(module_score, na.rm = TRUE)
  )

print(summary_df)

write.csv(
  summary_df,
  file = "AxonOutgrowth_ModuleScore_TP_vs_Combination_summary.csv",
  row.names = FALSE
)

############################################################
# Presynaptic terminal / vesicle release-related genes
# DotPlot + ModuleScore + boxplot + statistics
#
# 対象:
#   neuron_sub
#
# 色:
#   Combination = #2E8B57
#   TP          = #D55E00
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

# DotPlotで TP を上、Combination を下にしたいので
# factor順は Combination -> TP
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

# ==========================================================
# 1. geneset定義
# ==========================================================
presynaptic_genes <- c(
  "NRXN1", "PCLO", "BSN", "SNAP25", "RAB3A",
  "SYP", "SYN1", "STX1A", "VAMP2", "SV2A"
)

# 実際に存在する遺伝子のみ使う
presynaptic_genes_present <- presynaptic_genes[
  presynaptic_genes %in% rownames(neuron_sub)
]

cat("Genes found in neuron_sub:\n")
print(presynaptic_genes_present)

cat("Genes NOT found in neuron_sub:\n")
print(setdiff(presynaptic_genes, presynaptic_genes_present))

# ==========================================================
# 2. DotPlot
# ----------------------------------------------------------
# 右側凡例:
#   上 = Percent Expression
#   下 = Average Expression
# ==========================================================
p_dot_presynaptic <- DotPlot(
  object = neuron_sub,
  features = presynaptic_genes_present,
  group.by = "group",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Presynaptic terminal / vesicle release-related genes") +
  xlab("") +
  ylab("Group") +
  guides(
    size  = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_presynaptic)

ggsave(
  filename = "DotPlot_Presynaptic_VesicleRelease_genes_by_group.pdf",
  plot = p_dot_presynaptic,
  width = 9,
  height = 4.5
)

ggsave(
  filename = "DotPlot_Presynaptic_VesicleRelease_genes_by_group.png",
  plot = p_dot_presynaptic,
  width = 9,
  height = 4.5,
  dpi = 300
)

# ==========================================================
# 3. ModuleScore を計算
# ==========================================================
neuron_sub <- AddModuleScore(
  object = neuron_sub,
  features = list(presynaptic_genes_present),
  name = "PresynapticScore"
)

# 列名確認
cat("Metadata columns containing PresynapticScore:\n")
print(grep("PresynapticScore", colnames(neuron_sub@meta.data), value = TRUE))

score_col <- grep("PresynapticScore", colnames(neuron_sub@meta.data), value = TRUE)[1]

# ==========================================================
# 4. scoreの分布確認
# ==========================================================
VlnPlot(
  object = neuron_sub,
  features = score_col,
  group.by = "group",
  pt.size = 0.1
)

# ==========================================================
# 5. boxplot用データフレーム
# ==========================================================
plot_df_score <- neuron_sub@meta.data %>%
  dplyr::select(group, all_of(score_col)) %>%
  dplyr::rename(module_score = all_of(score_col))

head(plot_df_score)

# ==========================================================
# 6. 統計
# ==========================================================
wilcox_res <- wilcox.test(
  module_score ~ group,
  data = plot_df_score
)

ttest_res <- t.test(
  module_score ~ group,
  data = plot_df_score
)

cat("=== Wilcoxon test ===\n")
print(wilcox_res)

cat("=== t-test ===\n")
print(ttest_res)

stats_df <- data.frame(
  test = c("Wilcoxon rank-sum", "Welch t-test"),
  p_value = c(wilcox_res$p.value, ttest_res$p.value)
)

print(stats_df)

write.csv(
  stats_df,
  file = "Presynaptic_ModuleScore_TP_vs_Combination_statistics.csv",
  row.names = FALSE
)

# ==========================================================
# 7. boxplot
# ----------------------------------------------------------
# 図中に統計値は表示しない
# ==========================================================
plot_df_score$group <- factor(
  plot_df_score$group,
  levels = c("TP", "Combination")
)

group_colors <- c(
  "Combination" = "#2E8B57",
  "TP" = "#D55E00"
)

p_box_score <- ggplot(
  plot_df_score,
  aes(x = group, y = module_score, fill = group)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.18,
    size = 1.2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Presynaptic terminal / vesicle release-related module score",
    x = "",
    y = "Module score"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_box_score)

ggsave(
  filename = "Boxplot_Presynaptic_ModuleScore_TP_vs_Combination.pdf",
  plot = p_box_score,
  width = 6,
  height = 5
)

ggsave(
  filename = "Boxplot_Presynaptic_ModuleScore_TP_vs_Combination.png",
  plot = p_box_score,
  width = 6,
  height = 5,
  dpi = 300
)

# ==========================================================
# 8. 平均・中央値も保存
# ==========================================================
summary_df <- plot_df_score %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(module_score, na.rm = TRUE),
    sd = sd(module_score, na.rm = TRUE),
    median = median(module_score, na.rm = TRUE),
    IQR = IQR(module_score, na.rm = TRUE)
  )

print(summary_df)

write.csv(
  summary_df,
  file = "Presynaptic_ModuleScore_TP_vs_Combination_summary.csv",
  row.names = FALSE
)

############################################################
# Neurotransmission-related genes
# DotPlot + ModuleScore + boxplot + statistics
#
# 対象:
#   neuron_sub
#
# 色:
#   Combination = #2E8B57
#   TP          = #D55E00
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

# DotPlotで TP を上、Combination を下にしたいので
# factor順は Combination -> TP
neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("Combination", "TP")
)

# ==========================================================
# 1. geneset定義
# ==========================================================
neurotransmission_genes <- c(
  "TPH2", "DDC", "SLC18A2", "SLC17A8", "SLC22A3", "SLC6A4",
  "COMT", "MAOA", "MAOB", "ARRB2", "ADCYAP1", "TOR1A"
)

# 実際に存在する遺伝子のみ使う
neurotransmission_genes_present <- neurotransmission_genes[
  neurotransmission_genes %in% rownames(neuron_sub)
]

cat("Genes found in neuron_sub:\n")
print(neurotransmission_genes_present)

cat("Genes NOT found in neuron_sub:\n")
print(setdiff(neurotransmission_genes, neurotransmission_genes_present))

# ==========================================================
# 2. DotPlot
# ----------------------------------------------------------
# 右側凡例:
#   上 = Percent Expression
#   下 = Average Expression
# ==========================================================
p_dot_neurotransmission <- DotPlot(
  object = neuron_sub,
  features = neurotransmission_genes_present,
  group.by = "group",
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  ggtitle("Neurotransmission-related genes") +
  xlab("") +
  ylab("Group") +
  guides(
    size  = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_dot_neurotransmission)

ggsave(
  filename = "DotPlot_Neurotransmission_genes_by_group.pdf",
  plot = p_dot_neurotransmission,
  width = 10,
  height = 4.5
)

ggsave(
  filename = "DotPlot_Neurotransmission_genes_by_group.png",
  plot = p_dot_neurotransmission,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ==========================================================
# 3. ModuleScore を計算
# ==========================================================
neuron_sub <- AddModuleScore(
  object = neuron_sub,
  features = list(neurotransmission_genes_present),
  name = "NeurotransmissionScore"
)

# 列名確認
cat("Metadata columns containing NeurotransmissionScore:\n")
print(grep("NeurotransmissionScore", colnames(neuron_sub@meta.data), value = TRUE))

score_col <- grep("NeurotransmissionScore", colnames(neuron_sub@meta.data), value = TRUE)[1]

# ==========================================================
# 4. scoreの分布確認
# ==========================================================
VlnPlot(
  object = neuron_sub,
  features = score_col,
  group.by = "group",
  pt.size = 0.1
)

# ==========================================================
# 5. boxplot用データフレーム
# ==========================================================
plot_df_score <- neuron_sub@meta.data %>%
  dplyr::select(group, all_of(score_col)) %>%
  dplyr::rename(module_score = all_of(score_col))

head(plot_df_score)

# ==========================================================
# 6. 統計
# ==========================================================
wilcox_res <- wilcox.test(
  module_score ~ group,
  data = plot_df_score
)

ttest_res <- t.test(
  module_score ~ group,
  data = plot_df_score
)

cat("=== Wilcoxon test ===\n")
print(wilcox_res)

cat("=== t-test ===\n")
print(ttest_res)

stats_df <- data.frame(
  test = c("Wilcoxon rank-sum", "Welch t-test"),
  p_value = c(wilcox_res$p.value, ttest_res$p.value)
)

print(stats_df)

write.csv(
  stats_df,
  file = "Neurotransmission_ModuleScore_TP_vs_Combination_statistics.csv",
  row.names = FALSE
)

# ==========================================================
# 7. boxplot
# ----------------------------------------------------------
# 図中に統計値は表示しない
# ==========================================================
plot_df_score$group <- factor(
  plot_df_score$group,
  levels = c("TP", "Combination")
)

group_colors <- c(
  "Combination" = "#2E8B57",
  "TP" = "#D55E00"
)

p_box_score <- ggplot(
  plot_df_score,
  aes(x = group, y = module_score, fill = group)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.18,
    size = 1.2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Neurotransmission-related module score",
    x = "",
    y = "Module score"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_box_score)

ggsave(
  filename = "Boxplot_Neurotransmission_ModuleScore_TP_vs_Combination.pdf",
  plot = p_box_score,
  width = 6,
  height = 5
)

ggsave(
  filename = "Boxplot_Neurotransmission_ModuleScore_TP_vs_Combination.png",
  plot = p_box_score,
  width = 6,
  height = 5,
  dpi = 300
)

# ==========================================================
# 8. 平均・中央値も保存
# ==========================================================
summary_df <- plot_df_score %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(module_score, na.rm = TRUE),
    sd = sd(module_score, na.rm = TRUE),
    median = median(module_score, na.rm = TRUE),
    IQR = IQR(module_score, na.rm = TRUE)
  )

print(summary_df)

write.csv(
  summary_df,
  file = "Neurotransmission_ModuleScore_TP_vs_Combination_summary.csv",
  row.names = FALSE
)



############################################################
# 4 genesets:
# subclusterごとの ModuleScore boxplot
# x軸 = neuronal subcluster
# 各subcluster内で 左 = TP, 右 = Combination
# AddModuleScore列名を安全に取得する修正版
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ==========================================================
# 0. 前準備
# ==========================================================
DefaultAssay(neuron_sub) <- "RNA"

neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("TP", "Combination")
)

if (!"subcluster_name" %in% colnames(neuron_sub@meta.data)) {
  stop("neuron_sub@meta.data に subcluster_name がありません")
}

if (!is.factor(neuron_sub$subcluster_name)) {
  neuron_sub$subcluster_name <- factor(neuron_sub$subcluster_name)
}

group_colors <- c(
  "TP" = "#D55E00",
  "Combination" = "#2E8B57"
)

# ==========================================================
# 1. geneset定義
# ==========================================================
geneset_list <- list(
  AxonGuidance = c(
    "DCC", "NRP1", "NRP2", "PLXNA1",
    "ROBO1", "ROBO2", "EPHA4", "NRCAM",
    "FAT3", "IGF1R", "PTN", "RTN4RL1"
  ),
  AxonOutgrowth = c(
    "GAP43", "STMN2", "MAP1B", "DCX", "TUBA1A", "TUBB3"
  ),
  Presynaptic = c(
    "NRXN1", "PCLO", "BSN", "SNAP25", "RAB3A",
    "SYP", "SYN1", "STX1A", "VAMP2", "SV2A"
  ),
  Neurotransmission = c(
    "TPH2", "DDC", "SLC18A2", "SLC17A8", "SLC22A3", "SLC6A4",
    "COMT", "MAOA", "MAOB", "ARRB2", "ADCYAP1", "TOR1A"
  )
)

geneset_titles <- c(
  AxonGuidance = "Axon guidance / neurite pathfinding",
  AxonOutgrowth = "Axon outgrowth / regeneration-associated",
  Presynaptic = "Presynaptic terminal / vesicle release-related",
  Neurotransmission = "Neurotransmission-related"
)

# ==========================================================
# 2. AddModuleScoreを安全に実行し、追加列名を記録
# ==========================================================
score_name_map <- c()

for (set_name in names(geneset_list)) {
  
  gene_vec <- geneset_list[[set_name]]
  gene_present <- gene_vec[gene_vec %in% rownames(neuron_sub)]
  gene_missing <- setdiff(gene_vec, gene_present)
  
  cat("\n====================================\n")
  cat("Geneset:", set_name, "\n")
  cat("Present genes:\n")
  print(gene_present)
  cat("Missing genes:\n")
  print(gene_missing)
  
  if (length(gene_present) == 0) {
    warning(paste("No genes found for", set_name))
    score_name_map[set_name] <- NA_character_
    next
  }
  
  meta_before <- colnames(neuron_sub@meta.data)
  
  neuron_sub <- AddModuleScore(
    object = neuron_sub,
    features = list(gene_present),
    name = paste0(set_name, "Score")
  )
  
  meta_after <- colnames(neuron_sub@meta.data)
  new_cols <- setdiff(meta_after, meta_before)
  
  matched_new_cols <- new_cols[grepl(paste0("^", set_name, "Score"), new_cols)]
  
  if (length(matched_new_cols) == 0) {
    matched_all <- grep(
      paste0("^", set_name, "Score"),
      meta_after,
      value = TRUE
    )
    score_name_map[set_name] <- tail(matched_all, 1)
  } else {
    score_name_map[set_name] <- matched_new_cols[1]
  }
}

cat("\n=== score_name_map ===\n")
print(score_name_map)

# ==========================================================
# 3. 作図関数
# ==========================================================
make_subcluster_boxplot <- function(
    obj,
    score_col,
    plot_title,
    out_prefix,
    group_colors
) {
  
  if (is.na(score_col) || !score_col %in% colnames(obj@meta.data)) {
    stop(paste("score_col が存在しません:", score_col))
  }
  
  plot_df <- obj@meta.data %>%
    dplyr::select(subcluster_name, group, all_of(score_col)) %>%
    dplyr::rename(module_score = all_of(score_col)) %>%
    dplyr::filter(!is.na(module_score))
  
  plot_df$group <- factor(
    plot_df$group,
    levels = c("TP", "Combination")
  )
  
  summary_df <- plot_df %>%
    group_by(subcluster_name, group) %>%
    summarise(
      n = n(),
      mean = mean(module_score, na.rm = TRUE),
      sd = sd(module_score, na.rm = TRUE),
      median = median(module_score, na.rm = TRUE),
      IQR = IQR(module_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(
    summary_df,
    file = paste0(out_prefix, "_summary_by_subcluster.csv"),
    row.names = FALSE
  )
  
  stat_df <- plot_df %>%
    group_by(subcluster_name) %>%
    summarise(
      p_value_wilcox = tryCatch(
        wilcox.test(module_score ~ group)$p.value,
        error = function(e) NA_real_
      ),
      p_value_ttest = tryCatch(
        t.test(module_score ~ group)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  write.csv(
    stat_df,
    file = paste0(out_prefix, "_statistics_by_subcluster.csv"),
    row.names = FALSE
  )
  
  p <- ggplot(
    plot_df,
    aes(x = subcluster_name, y = module_score, fill = group)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.85
    ) +
    geom_jitter(
      aes(color = group),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.75
      ),
      size = 0.8,
      alpha = 0.45,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(
      title = plot_title,
      x = "Neuronal subcluster",
      y = "Module score"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  print(p)
  
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p,
    width = 10,
    height = 5
  )
  
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  return(list(
    plot = p,
    summary = summary_df,
    stats = stat_df
  ))
}

# ==========================================================
# 4. 4 geneset それぞれで実行
# ==========================================================
res_AxonGuidance <- make_subcluster_boxplot(
  obj = neuron_sub,
  score_col = unname(score_name_map["AxonGuidance"]),
  plot_title = geneset_titles["AxonGuidance"],
  out_prefix = "SubclusterBox_AxonGuidance_ModuleScore",
  group_colors = group_colors
)

res_AxonOutgrowth <- make_subcluster_boxplot(
  obj = neuron_sub,
  score_col = unname(score_name_map["AxonOutgrowth"]),
  plot_title = geneset_titles["AxonOutgrowth"],
  out_prefix = "SubclusterBox_AxonOutgrowth_ModuleScore",
  group_colors = group_colors
)

res_Presynaptic <- make_subcluster_boxplot(
  obj = neuron_sub,
  score_col = unname(score_name_map["Presynaptic"]),
  plot_title = geneset_titles["Presynaptic"],
  out_prefix = "SubclusterBox_Presynaptic_ModuleScore",
  group_colors = group_colors
)

res_Neurotransmission <- make_subcluster_boxplot(
  obj = neuron_sub,
  score_col = unname(score_name_map["Neurotransmission"]),
  plot_title = geneset_titles["Neurotransmission"],
  out_prefix = "SubclusterBox_Neurotransmission_ModuleScore",
  group_colors = group_colors
)


############################################################
# 各 subcluster ごとに
# Axon Guidance / Axon Outgrowth / Presynaptic /
# Neurotransmission score の TP vs Combination 統計
############################################################

library(dplyr)

# ==========================================================
# 0. 前提確認
# ==========================================================
if (!"subcluster_name" %in% colnames(neuron_sub@meta.data)) {
  stop("neuron_sub@meta.data に subcluster_name がありません")
}

if (!"group" %in% colnames(neuron_sub@meta.data)) {
  stop("neuron_sub@meta.data に group がありません")
}

neuron_sub$group <- factor(
  neuron_sub$group,
  levels = c("TP", "Combination")
)

# ==========================================================
# 1. score列名を自動取得
# ----------------------------------------------------------
# AddModuleScoreを複数回回していても、最後の列を拾う
# ==========================================================
score_name_map <- c(
  AxonGuidance = tail(grep("^AxonGuidanceScore", colnames(neuron_sub@meta.data), value = TRUE), 1),
  AxonOutgrowth = tail(grep("^AxonOutgrowthScore", colnames(neuron_sub@meta.data), value = TRUE), 1),
  Presynaptic = tail(grep("^PresynapticScore", colnames(neuron_sub@meta.data), value = TRUE), 1),
  Neurotransmission = tail(grep("^NeurotransmissionScore", colnames(neuron_sub@meta.data), value = TRUE), 1)
)

cat("Detected score columns:\n")
print(score_name_map)

# ==========================================================
# 2. 各 score について subclusterごとの統計を出す関数
# ==========================================================
calc_subcluster_stats <- function(obj, score_col, score_label, out_prefix) {
  
  if (is.na(score_col) || !score_col %in% colnames(obj@meta.data)) {
    stop(paste("score_col が存在しません:", score_col))
  }
  
  plot_df <- obj@meta.data %>%
    dplyr::select(subcluster_name, group, all_of(score_col)) %>%
    dplyr::rename(module_score = all_of(score_col)) %>%
    dplyr::filter(!is.na(module_score))
  
  plot_df$group <- factor(
    plot_df$group,
    levels = c("TP", "Combination")
  )
  
  stat_df <- plot_df %>%
    dplyr::group_by(subcluster_name) %>%
    dplyr::summarise(
      n_TP = sum(group == "TP"),
      n_Combination = sum(group == "Combination"),
      p_value_wilcox = tryCatch(
        wilcox.test(module_score ~ group)$p.value,
        error = function(e) NA_real_
      ),
      p_value_ttest = tryCatch(
        t.test(module_score ~ group)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(p_value_wilcox)
  
  cat("\n====================================\n")
  cat(score_label, "\n")
  print(stat_df)
  
  write.csv(
    stat_df,
    file = paste0(out_prefix, "_statistics_by_subcluster.csv"),
    row.names = FALSE
  )
  
  return(stat_df)
}

# ==========================================================
# 3. 4種類の score で実行
# ==========================================================
stat_axon_guidance <- calc_subcluster_stats(
  obj = neuron_sub,
  score_col = unname(score_name_map["AxonGuidance"]),
  score_label = "Axon Guidance score",
  out_prefix = "Subcluster_AxonGuidance"
)

stat_axon_outgrowth <- calc_subcluster_stats(
  obj = neuron_sub,
  score_col = unname(score_name_map["AxonOutgrowth"]),
  score_label = "Axon Outgrowth score",
  out_prefix = "Subcluster_AxonOutgrowth"
)

stat_presynaptic <- calc_subcluster_stats(
  obj = neuron_sub,
  score_col = unname(score_name_map["Presynaptic"]),
  score_label = "Presynaptic score",
  out_prefix = "Subcluster_Presynaptic"
)

stat_neurotransmission <- calc_subcluster_stats(
  obj = neuron_sub,
  score_col = unname(score_name_map["Neurotransmission"]),
  score_label = "Neurotransmission score",
  out_prefix = "Subcluster_Neurotransmission"
)

# ==========================================================
# 4. 見やすく1つにまとめる
# ==========================================================
stat_axon_guidance$score_type <- "Axon Guidance"
stat_axon_outgrowth$score_type <- "Axon Outgrowth"
stat_presynaptic$score_type <- "Presynaptic"
stat_neurotransmission$score_type <- "Neurotransmission"

stat_all <- dplyr::bind_rows(
  stat_axon_guidance,
  stat_axon_outgrowth,
  stat_presynaptic,
  stat_neurotransmission
) %>%
  dplyr::select(score_type, everything())

cat("\n====================================\n")
cat("All score statistics combined:\n")
print(stat_all)

write.csv(
  stat_all,
  file = "Subcluster_AllScores_statistics_by_subcluster.csv",
  row.names = FALSE
)