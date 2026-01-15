library(Matrix);library(Seurat);library(data.table)
rm(list = ls())

setwd()
#####  1.load  #################################################################
####
# GSE127465
raw_dir <- "/GSE127465_RAW/"
meta_file <- "/GSE127465_RAW/GSE127465_mouse_cell_metadata_15939x12.tsv.gz"
Meta_all <- fread(meta_file, data.table = FALSE)
raw_files <- list.files(raw_dir,pattern = "_raw_counts.tsv.gz$",full.names = TRUE)
make_seurat_from_raw <- function(file, Meta_all) {
  message("Processing: ", basename(file))
  expr <- fread(file, data.table = FALSE)
  barcodes <- expr$barcode
  mat <- as.matrix(expr[, -1])
  mat <- t(mat)  # gene x cell
  rownames(mat) <- colnames(expr)[-1]
  sample_id <- sub("_raw_counts.tsv.gz", "", basename(file))
  lib_id <- sub(".*_mouse_", "", sample_id)
  cell_ids <- paste(lib_id, barcodes, sep = "_")
  colnames(mat) <- cell_ids
  mat <- Matrix(mat, sparse = TRUE)
  seu <- CreateSeuratObject(counts = mat,project = sample_id)
  Meta_sub <- Meta_all[Meta_all$Library == lib_id, ]
  Meta_sub$cell_id <- paste(Meta_sub$Library, Meta_sub$Barcode, sep = "_")
  common_cells <- intersect(colnames(seu), Meta_sub$cell_id)
  seu <- subset(seu, cells = common_cells)
  stopifnot(all(Meta_sub$cell_id == colnames(seu)))
  rownames(Meta_sub) <- Meta_sub$cell_id
  seu <- AddMetaData(seu, Meta_sub)
  return(seu)
}
seu_list <- lapply(raw_files, make_seurat_from_raw, Meta_all = Meta_all)
names(seu_list) <- sub("_raw_counts.tsv.gz", "", basename(raw_files))
GSE127465 <- merge(x = seu_list[[1]],y = seu_list[-1],add.cell.ids = names(seu_list),project = "GSE127465")
GSE127465$Samples <- GSE127465$Library
GSE127465$Source <- "GSE127465"
GSE127465$orig.ident <- GSE127465$cell_id

# GSE133604
raw_dir <- "/GSE133604"
genes <- fread(file.path(raw_dir, "GSE133604_genes.tsv.gz"),header = FALSE)[,2]
mtx_files <- list.files(raw_dir,pattern = "_matrix.mtx.gz$",full.names = TRUE)
make_seurat_from_shared_mtx <- function(mtx_file, genes, barcodes) {
  message("Processing: ", basename(mtx_file))
  mat <- readMM(mtx_file)
  dim(mat)
  sample_id <- sub("_matrix.mtx.gz", "", basename(mtx_file))
  rownames(mat) <- make.unique(genes$V2)
  colnames(mat) <- paste0(sample_id, "_cell", seq_len(ncol(mat)))
  seu <- CreateSeuratObject(counts = mat,project = "GSE133604")
  seu$sample_id <- sample_id
  return(seu)
}
seu_list_133604 <- lapply(mtx_files,make_seurat_from_shared_mtx,genes = genes)
GSE133604 <-  merge(seu_list_133604[[1]],y = seu_list_133604[-1])
GSE133604$Samples <- GSE133604$sample_id
GSE133604$Source <- "GSE133604"
GSE133604$orig.ident <- rownames(GSE133604@meta.data)

# GSE134255
make_seurat_from_counts_csv <- function(file, project = "GSE134255") {
  message("Processing: ", basename(file))
  expr <- fread(file, data.table = FALSE)
  gene_names <- expr[[1]]
  counts <- as.matrix(expr[, -1])
  rownames(counts) <- gene_names
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  sample_id <- sub("_counts.csv.gz$", "", basename(file))
  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,project = project)
  seu$Sample <- sample_id
  return(seu)
}
raw_dir <- "/GSE134255/"
mtx_files <- list.files(raw_dir,pattern = "_counts.csv.gz$",full.names = TRUE)
seu_list_134255 <- lapply(mtx_files,make_seurat_from_counts_csv)
GSE134255 <-  merge(seu_list_134255[[1]],y = seu_list_134255[-1])
GSE134255$Samples <- GSE134255$Sample
GSE134255$Source <- "GSE134255"
GSE134255$orig.ident <- rownames(GSE134255@meta.data)
# GSE136206
make_seurat_from_10x_nonstandard <- function(mtx_file, project = "GSE136206") {
  prefix <- sub("\\.matrix\\.mtx\\.gz$", "", basename(mtx_file))
  dir <- dirname(mtx_file)
  message("Processing: ", prefix)
  feat_file <- file.path(dir, paste0(prefix, ".features.tsv.gz"))
  bc_file   <- file.path(dir, paste0(prefix, ".barcodes.tsv.gz"))
  stopifnot(file.exists(mtx_file))
  stopifnot(file.exists(feat_file))
  stopifnot(file.exists(bc_file))
  mat <- ReadMtx(mtx = mtx_file,features = feat_file,cells = bc_file,feature.column = 2)
  colnames(mat) <- paste(prefix, colnames(mat), sep = "_")
  seu <- CreateSeuratObject(counts = mat,project = project)
  seu$Sample <- prefix
  return(seu)
}
raw_dir <- "/GSE136206"
mtx_files <- list.files(raw_dir,pattern = "matrix.mtx.gz$",full.names = TRUE)
seu_list_136206 <- lapply(mtx_files,make_seurat_from_10x_nonstandard)
names(seu_list_136206) <- sub("\\.matrix\\.mtx\\.gz$", "", basename(mtx_files))
GSE136206 <-  merge(seu_list_136206[[1]],y = seu_list_136206[-1])
GSE136206$Samples <- GSE136206$Sample
GSE136206$Source <- "GSE136206"
GSE136206$orig.ident <- rownames(GSE136206@meta.data)
# GSE157561
h5_files <- list.files(path = "/GSE157561",pattern = "\\.h5$",full.names = TRUE)
seu_list_157561 <- lapply(h5_files, function(f) {
  prefix <- sub("\\.h5$", "", basename(f))
  counts <- Read10X_h5(f)
  colnames(counts) <- paste(prefix, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts, project = "GSE157561")
  seu$Sample <- prefix
  return(seu)
})
names(seu_list_157561) <- sub("\\.h5$", "", basename(h5_files))
GSE157561 <-  merge(seu_list_157561[[1]],y = seu_list_157561[-1])
GSE157561$Samples <- GSE157561$Sample
GSE157561$Source <- "GSE157561"
GSE157561$orig.ident <- rownames(GSE157561@meta.data)
GSE157561@meta.data <- GSE157561@meta.data[,metadata_col]

metadata_col <- c("orig.ident","nCount_RNA","nFeature_RNA","Samples","Source")
GSE127465@meta.data <- GSE127465@meta.data[,metadata_col]
GSE133604@meta.data <- GSE133604@meta.data[,metadata_col]
GSE134255@meta.data <- GSE134255@meta.data[,metadata_col]
GSE136206@meta.data <- GSE136206@meta.data[,metadata_col]
GSE157561@meta.data <- GSE157561@meta.data[,metadata_col]
## Merge
Tumor_AllCell <- merge(x = GSE127465,y =list(GSE133604,GSE134255,GSE136206,GSE157561), merge.data=T,project = "SeuratProject")


counts <- GetAssayData(Tumor_AllCell, slot = "counts")
keep_genes <- rowSums(counts > 0) >= 10
Tumor_AllCell <- Tumor_AllCell[keep_genes, ]
Tumor_AllCell[["percent.mt"]] <- PercentageFeatureSet(Tumor_AllCell,pattern = "^mt-")
Tumor_AllCell <- subset(Tumor_AllCell,subset =nCount_RNA >= 600 & nCount_RNA <= 120000 &
                          nFeature_RNA >= 400 & nFeature_RNA <= 8000 &percent.mt <= 10)



obj_M <- NormalizeData(Tumor_AllCell, normalization.method = "LogNormalize", scale.factor = 1e4,verbose = FALSE)
obj_M <- FindVariableFeatures(obj_M, selection.method = "vst",nfeatures = 3000, verbose = FALSE)
obj_M <- ScaleData(object = obj_M, verbose = FALSE)
obj_M <- RunPCA(obj_M, npcs = 30, verbose = FALSE)
obj_M <- RunHarmony(obj_M,group.by.vars = "Source",dims.use = 1:30)
obj_M <- FindNeighbors(obj_M, reduction = "harmony",dims = 1:30, verbose = FALSE)
obj_M <- FindVariableFeatures(obj_M, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
obj_M <- FindClusters(obj_M, resolution = 0.8,algorithm = 1, verbose = FALSE)

ref <- readRDS("/Mouse_ref.rds")
sce <- as.SingleCellExperiment(obj_M)
clusters <- obj_M$seurat_clusters
stopifnot(length(clusters) == ncol(sce))

pred <- SingleR(test = sce,ref = ref,labels = ref$label.main,clusters = clusters,
                assay.type.test = "logcounts",assay.type.ref  = "logcounts")

cluster_anno <- data.frame(cluster = rownames(pred),label   = pred$labels,stringsAsFactors = FALSE)
cluster_anno$label_pruned <- pred$pruned.labels
stopifnot(length(clusters) == ncol(sce))

m <- setNames(pred$labels, rownames(pred)) 
sce$SingleR_label_main <- unname(m[as.character(clusters)])
m2 <- setNames(pred$pruned.labels, rownames(pred))
sce$SingleR_label_main_pruned <- unname(m2[as.character(clusters)])
obj_M$SingleR_label_main <- sce$SingleR_label_main[match(colnames(obj_M), colnames(sce))]

Colors <- c("#daa4c0","#e99c8d","#d19959","#7789b1","#e4b769","#6596ae","#649d59","#a26796","#9f8baf","#b62f2e","#d0717f")

Idents(obj_M) <- obj_M$Maintypes
p <- DimPlot(obj_M,reduction = "umap",group.by = "Maintypes",pt.size = 1,cols = Colors)
p$layers[[1]]$aes_params$alpha <- 0.3
p +labs(title = "Inhouse cohort Landscape") +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.title = element_blank(),legend.text = element_text(size = 9),legend.position = "bottom") +
  coord_fixed()



genes_of_interest <- c("Comt", "Maoa", "Maob")
expr_mat <- FetchData(obj_M, vars = c(genes_of_interest, "Maintypes"))
expr_long <- expr_mat %>%pivot_longer(cols = all_of(genes_of_interest),names_to = "Gene",values_to = "Expression")

avg_expr <- expr_long %>%group_by(Maintypes, Gene) %>%summarise(Expression = mean(Expression), .groups = "drop")
cluster_order <- arrange(avg_expr[avg_expr$Gene=="COMT",],Expression)

ggplot(avg_expr, aes(x = Gene, y = Maintypes, fill = Expression)) +geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(avg_expr$Expression)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank()) +labs(fill = "Avg Expression")+
  scale_y_discrete(limits = cluster_order$Maintypes) 




