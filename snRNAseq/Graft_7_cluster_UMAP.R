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