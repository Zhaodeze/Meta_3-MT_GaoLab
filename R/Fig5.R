library(Seurat);library(SingleR);library(SingleCellExperiment);library(celldex);library(harmony);
rm(list = ls())

#####  1.load  #################################################################
####
metadata_col <- c("orig.ident","nCount_RNA","nFeature_RNA","Samples","TissueType","Source")
## CRA001160
CRA001160 <- readRDS('/CRA001160_WholeTissue.rds.gz')
CRA001160@meta.data <-  CRA001160@meta.data%>%mutate(TissueType=ifelse(grepl("^N", Samples), "Normal",ifelse(grepl("^T", Samples), "Tumor", NA)))
CRA001160$Source <- "CRA001160"
CRA001160$orig.ident <- colnames(CRA001160)
CRA001160@meta.data <- CRA001160@meta.data[,metadata_col]
CRA001160_T <- subset(CRA001160,TissueType=="Tumor")
## GSE123813
GSE123813 <- readRDS('/GSE123813_bcc_wholetissue.rds.gz')
colnames(GSE123813@meta.data)[5] <- "Samples"
GSE123813$TissueType <- "Tumor"
GSE123813$Source <- "GSE123813"
GSE123813$orig.ident <- colnames(GSE123813)
GSE123813@meta.data <- GSE123813@meta.data[,metadata_col]
GSE123813_T <- subset(GSE123813,TissueType=="Tumor")
## GSE149614
GSE149614 <- readRDS('/GSE149614_HCC_wholetissue.rds.gz')
colnames(GSE149614@meta.data)
colnames(GSE149614@meta.data)[5] <- "Samples"
colnames(GSE149614@meta.data)[7] <- "TissueType"
GSE149614$Source <- "GSE149614"
GSE149614$orig.ident <- colnames(GSE149614)
GSE149614@meta.data <- GSE149614@meta.data[,metadata_col]
GSE149614_T <- subset(GSE149614,TissueType=="Tumor")
## GSE176078
GSE176078 <- readRDS('/GSE176078_WholeTissue.rds.gz')
GSE176078$Samples <- GSE176078$orig.ident
GSE176078$TissueType <- "Tumor"
GSE176078$Source <- "GSE176078"
GSE176078$orig.ident <- colnames(GSE176078)
GSE176078@meta.data <- GSE176078@meta.data[,metadata_col]
GSE176078_T <- subset(GSE176078,TissueType=="Tumor")
## GSE185344
GSE185344 <- readRDS('/GSE185344_PH_scRNA.final.rds')
GSE185344 <- GSE185344$obj
GSE185344$Samples <- sub("^(HYW_[0-9]+)_.*", "\\1", GSE185344$orig.ident)
GSE185344$TissueType <- sub("^HYW_[0-9]+_([^_]+).*", "\\1", GSE185344$orig.ident)
GSE185344$Source <- "GSE185344"
GSE185344$orig.ident <- colnames(GSE185344)
GSE185344@meta.data <- GSE185344@meta.data[,metadata_col]
GSE185344_T <- subset(GSE185344,TissueType=="Tumor")
## GSE205506
GSE205506 <- readRDS('/GSE205506_WholeTissue.rds.gz')
GSE205506_info <- data.table::fread('/GSE205506_Info.csv')%>%
  select(GEO_ID, Samples, TissueType) %>%distinct()
GSE205506$Samples <- GSE205506_info$Samples[match(GSE205506$orig.ident, GSE205506_info$GEO_ID)]
GSE205506$TissueType <- GSE205506_info$TissueType[match(GSE205506$orig.ident, GSE205506_info$GEO_ID)]
GSE205506$Source <- "GSE205506"
GSE205506$orig.ident <- colnames(GSE205506)
GSE205506@meta.data <- GSE205506@meta.data[,metadata_col]
GSE205506_T <- subset(GSE205506,TissueType=="Tumor")
## GSE215121
GSE215121 <- readRDS('/GSE215121_WholeTissue_object.rds.gz')
GSE215121$Samples <- sub("^GSM[0-9]+_", "", GSE215121$orig.ident)
GSE215121$TissueType <- "Tumor"
GSE215121$Source <- "GSE215121"
GSE215121$orig.ident <- colnames(GSE215121)
GSE215121@meta.data <- GSE215121@meta.data[,metadata_col]
GSE215121_T <- subset(GSE215121,TissueType=="Tumor")
## EMTAB6149
EMTAB6149 <- readRDS('/EMTAB6149_NSCLC_seurat.rds.gz')
EMTAB6149$Samples <- EMTAB6149$orig.ident
colnames(EMTAB6149@meta.data)[13] <- "TissueType"
EMTAB6149$Source <- "EMTAB6149"
EMTAB6149$orig.ident <- paste(EMTAB6149$Samples, colnames(EMTAB6149), sep = "_")
EMTAB6149@meta.data <- EMTAB6149@meta.data[,metadata_col]
EMTAB6149_T <- subset(EMTAB6149,TissueType=="Tumor")
## Merge
Tumor_AllCell <- merge(x = CRA001160_T, y =list(GSE123813_T,GSE149614_T,GSE176078_T,GSE185344_T,GSE205506_T,GSE215121_T,EMTAB6149_T), merge.data=T,project = "SeuratProject")


counts <- GetAssayData(Tumor_AllCell, slot = "counts")
keep_genes <- rowSums(counts > 0) >= 10
Tumor_AllCell <- Tumor_AllCell[keep_genes, ]
Tumor_AllCell[["percent.mt"]] <- PercentageFeatureSet(Tumor_AllCell,pattern = "^MT-")
Tumor_AllCell <- subset(Tumor_AllCell,subset =nCount_RNA >= 600 & nCount_RNA <= 120000 &
                        nFeature_RNA >= 400 & nFeature_RNA <= 8000 &percent.mt <= 10)


obj <- NormalizeData(Tumor_AllCell, normalization.method = "LogNormalize", scale.factor = 1e4,verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(object = obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunHarmony(obj,group.by.vars = "Source",dims.use = 1:30)
obj <- FindNeighbors(obj, reduction = "harmony",dims = 1:30, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8,algorithm = 1, verbose = FALSE)

ref <- readRDS("/HPCA_ref.rds")
sce <- as.SingleCellExperiment(obj)
clusters <- obj$seurat_clusters
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
obj$SingleR_label_main <- sce$SingleR_label_main[match(colnames(obj), colnames(sce))]

Colors <- c("#daa4c0","#e99c8d","#d19959","#7789b1","#e4b769","#6596ae","#649d59","#a26796","#9f8baf","#b62f2e","#d0717f")

Idents(obj) <- obj$Maintypes
p <- DimPlot(obj,reduction = "umap",group.by = "Maintypes",pt.size = 1,cols = Colors)
p$layers[[1]]$aes_params$alpha <- 0.3
p +labs(title = "Inhouse cohort Landscape") +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
    panel.grid = element_blank(),plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_blank(),legend.text = element_text(size = 9),legend.position = "bottom") +
  coord_fixed()



genes_of_interest <- c("COMT", "MAOA", "MAOB")
expr_mat <- FetchData(obj, vars = c(genes_of_interest, "Maintypes"))
expr_long <- expr_mat %>%pivot_longer(cols = all_of(genes_of_interest),names_to = "Gene",values_to = "Expression")

avg_expr <- expr_long %>%group_by(Maintypes, Gene) %>%summarise(Expression = mean(Expression), .groups = "drop")
cluster_order <- arrange(avg_expr[avg_expr$Gene=="COMT",],Expression)

ggplot(avg_expr, aes(x = Gene, y = Maintypes, fill = Expression)) +geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(avg_expr$Expression)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank()) +labs(fill = "Avg Expression")+
  scale_y_discrete(limits = cluster_order$Maintypes) 

