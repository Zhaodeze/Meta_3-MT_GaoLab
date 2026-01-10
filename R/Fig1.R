####################################
####  Fig1 Prepare #################
####################################
rm(list = ls())
library(ggplot2);library(dplyr);library(magrittr);library(tidyr);library(matrixStats);library(circlize);library(GSVA);#ver.1.40.1
library(RColorBrewer);library(ComplexHeatmap);library(survival);library(reshape2);library(survminer);library(NMF);
#Fig1A  Meta_TvsN_All
#Fig1B  Impact_Pathway
#Fig1C  sample.TN,p_hsa,Meta_109,Clinic_109,Data_GSVA,KEGG_Catogorize
#Fig1D  Cluster_Compare,WES_index_sp,Path_TvsN_All,WES_Order
#Fig1E  Meta_109_Raw
#Fig1F  Sample.info,Meta_Class

load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig1/Fig1_All_DataMatrix.Rdata")
####################################
####  Fig1A  #######################
####################################
colors<-  c('#81cac0','#b1b0d0','#f37d74','#77a8c6','#a6d177','#f3c2da','#ac7cb7','#f5a86a','grey')
####   Draw   ####
ggplot(Meta_TvsN_All, aes(x = Log2FC, y = -log10(FDR)))+geom_point(aes(color = Color_Category),size=5) +
  scale_color_manual(values = colors) +geom_vline(xintercept= 0,lty=4,col="Moccasin",lwd=2)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="Moccasin",lwd=2) +theme_bw(base_size = 16)+ 
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=20),axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  labs(x="log2 Fold Change between tumor and adjacent",y="-log10 (p-value)")+xlim(c(-max(abs(Meta_TvsN_All$Log2FC)),max(abs(Meta_TvsN_All$Log2FC))))

####################################
####  Fig1B  #######################
####################################
ggplot(Impact_Pathway,aes(x = Impact,y = log10P,color = group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_vline(xintercept = 0,linetype = 4,colour = "grey40",linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05),linetype = 4,colour = "grey40",linewidth = 0.8) +
  ggrepel::geom_text_repel(data = subset(Impact_Pathway, group == "sig"),aes(label = pathway_name),
                           size = 5,box.padding = 0.4,point.padding = 0.3,max.overlaps = Inf,show.legend = FALSE) +
  scale_color_manual(values = c("grey", "#FF3333")) +
  labs(x = "Impact",y = "Log10(p)",title = "Pathway Impact vs Significance") +
  theme_classic(base_size = 18)
####################################
####  Fig1C  #######################
####################################
#### NMF Data Prepare #######
Data_GSVA <- GSVA::gsva(as.matrix(Meta_109[,-1]),p_hsa,method = 'gsva',min.sz=5,mx.diff=TRUE,verbose=F,parallel.sz=8)%>%
  {rownames(.) <- gsub(" \\- Homo sapiens \\(human\\)", "", rownames(.));df <- t(.) %>% as.data.frame();df[] <- lapply(df, as.numeric);
  df$ID<- rownames(df);df$Type  <- sample.TN$group[match(df$ID, sample.TN$ID)];df$Order <- gsub("M_", "", df$ID) %>% as.numeric();df %>% dplyr::arrange(Order);df$Pair  <- rep(1:109, each = 2);df}
#### Calculate Diff Pathway ####
Index <- head(colnames(Data_GSVA), -4)
Path_TvsN_All <- purrr::map_dfr(Index, function(sub) { Data <- Data_GSVA %>%dplyr::select(ID, Type, Pair, FPKM = all_of(sub)) %>%
  filter(!is.na(FPKM)) %>%dplyr::arrange(Type, Pair);N_FPKM <- Data$FPKM[Data$Type == "Non_Tumor"];T_FPKM <- Data$FPKM[Data$Type == "Tumor"]
  tibble(Symbol = sub,Difference = mean(Data$FPKM[Data$Type == "Tumor"], na.rm = TRUE)-mean(Data$FPKM[Data$Type == "Non_Tumor"], na.rm = TRUE),Pvalue = wilcox.test(T_FPKM, N_FPKM, paired = TRUE)$p.value)})%>%
  mutate(FDR = p.adjust(.$Pvalue,method="fdr"),Catogorize= KEGG_Catogorize$Catogorize[match(Symbol,KEGG_Catogorize$Pathway)])%>%
  mutate(Catogorize = dplyr::case_when(Catogorize %in% c("Amino acid metabolism", "Metabolism of other amino acids") ~ "Amino acid",Catogorize == "Xenobiotics biodegradation and metabolism" ~ "Drug",Catogorize == "Carbohydrate metabolism" ~ "Carbohydrate",
                                       Catogorize == "Lipid metabolism" ~ "Lipid",TRUE ~ "Other")) %>%dplyr::arrange(.,FDR)
#### NMF ####
zero_pro <- function(x){x <- length(which(x == 0))/length(x);return(x)}
TME_ssGSESA_new <- Data_GSVA %>%{ set_rownames(., .$ID) } %>%.[sample.TN[group == "Tumor", ]$ID,filter(Path_TvsN_All, FDR < 0.05)$Symbol] %>%t() %>%
  {mat1 <- .;mat2 <- .;mat1[which(mat1 < 0)] <- 0;mat2[which(mat1 > 0)] <- 0;mat2 <- abs(mat2);rownames(mat2) <- paste0(rownames(mat2), ".1");mat_new <- rbind(mat1, mat2);zero_pro_re <- apply(mat_new, 1, zero_pro);mat_new[which(zero_pro_re < 1), ]}
seed= 2066;nmf_res <- nmf(x=TME_ssGSESA_new, rank=2:5, seed=seed, nrun = 100);coph <- nmf_res$measures$cophenetic;opt_rank_name <- as.character(which.max(coph[-length(coph)]-coph[-1]));nmf_res_opt <- nmf_res$fit[[opt_rank_name]];
feat_opt_idx <- extractFeatures(object=nmf_res_opt, method = "max", format = "list", nodups = TRUE);feat_opt <- rownames(TME_ssGSESA_new)[as.numeric(na.omit(unique(unlist(feat_opt_idx))))]
nmf_cluster_opt <- nmf(TME_ssGSESA_new[feat_opt, ], rank = as.numeric(opt_rank_name), seed = seed)
Clinic_109$Cluster  <- paste0('M',predict(nmf_cluster_opt))[match(Clinic_109$ID,colnames(TME_ssGSESA_new))]
#### KM plot ####
fit <- survival::survfit(survival::Surv( Overall_survial, Survial_status) ~  Cluster,
                         data=Clinic_109, na.action=na.exclude)
model <-        survdiff(survival::Surv( Overall_survial, Survial_status) ~  Cluster,
                         data=Clinic_109, na.action=na.exclude)
KMP <- 1-pchisq(model$chisq, df=length(levels(factor(Clinic_109$Cluster)))-1)
sur_p <- ggsurvplot(fit, data = Clinic_109,legend.title="NMF group",
                    pval = TRUE, risk.table = T,palette = c("#db5243","#4db1c7","#309d88"),
                    xlab="Months",ylab="Overall Survival rate (OS)")
colnames(Clinic_109)
####################################
####  Fig1D  #######################
####################################
#### Cluster_Compare ####
Cluster_Compare <- Clinic_109[,c( "T_ID","Cluster","Cluster_Gao",                                
                                  "Cluster_FuchuHe","Cluster_XiaolongLiu")]%>%
  mutate(Proteomic_sub   = dplyr::case_when(Cluster_Gao %in% c("SI") ~ "S-Mb",TRUE ~ "Other"),
         FuchuHe_sub   = dplyr::case_when(Cluster_FuchuHe %in% c("SI") ~ "SI",TRUE ~ "Other"),
         XiaolongLiu_sub   = dplyr::case_when(Cluster_XiaolongLiu %in% c("SI") ~ "SI",TRUE ~ "Other"))
Chisq <- Cluster_Compare[,c( "Cluster","Proteomic_sub",                                
                             "FuchuHe_sub","XiaolongLiu_sub")] %>%set_rownames(Cluster_Compare$T_ID)%>%t()%>%cbind(Features=rownames(.),.)%>%as.data.frame()
Result_Chisq <- cbind(Features=Chisq$Features,P=sapply(Chisq$Features,function(sub){sample.Tumor_sub <- Cluster_Compare[,c("Cluster",sub)]
colnames(sample.Tumor_sub)[2] <- "candidate";contingency_table <- table(sample.Tumor_sub$Cluster, sample.Tumor_sub$candidate);chisq.test(contingency_table)$p.value}))

#### Cluster_Compare M1 ####
Cluster_Compare_M1 <- Cluster_Compare[which(Cluster_Compare$Cluster=="M1"),-c(3:5)]
Cluster_Compare_M1$T_ID <- factor(Cluster_Compare_M1$T_ID,levels = WES_M1_Order)
Cluster_Compare_M1_melt <-  reshape2::melt(Cluster_Compare_M1, id.vars = c("T_ID"))

pM1 <- ggplot(Cluster_Compare_M1_melt, aes(x = T_ID, y = variable)) +geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("M1"="#db5243","M2"="#4db1c7","M3"="#309d88","S-Mb"="#9dc591","S-Pf"="lightgrey","S-Me"="lightgrey",
                               "SI"="#f9c794","SII"="lightgrey","SIII"="lightgrey","S1"="#9492c2","S2"="lightgrey","S3"="lightgrey")) +
  theme_minimal() +labs(x = "Samples", y = "Categories", fill = "Values") +theme(axis.text.x = element_text(angle = 45, hjust = 1))
#### Cluster_Compare M2 ####
Cluster_Compare_M2 <- Cluster_Compare[which(Cluster_Compare$Cluster=="M2"),-c(3:5)]
Cluster_Compare_M2$T_ID <- factor(Cluster_Compare_M2$T_ID,levels = WES_M2_Order)
Cluster_Compare_M2_melt <-  reshape2::melt(Cluster_Compare_M2, id.vars = c("T_ID"))

pM2 <- ggplot(Cluster_Compare_M2_melt, aes(x = T_ID, y = variable)) +geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("M1"="#db5243","M2"="#4db1c7","M3"="#309d88","S-Mb"="#9dc591","S-Pf"="lightgrey","S-Me"="lightgrey",
                               "SI"="#f9c794","SII"="lightgrey","SIII"="lightgrey","S1"="#9492c2","S2"="lightgrey","S3"="lightgrey")) +
  theme_minimal() +labs(x = "Samples", y = "Categories", fill = "Values") +theme(axis.text.x = element_text(angle = 45, hjust = 1))
#### Cluster_Compare M3 ####
Cluster_Compare_M3 <- Cluster_Compare[which(Cluster_Compare$Cluster=="M3"),-c(3:5)]
Cluster_Compare_M3$T_ID <- factor(Cluster_Compare_M3$T_ID,levels = WES_M3_Order)
Cluster_Compare_M3_melt <-  reshape2::melt(Cluster_Compare_M3, id.vars = c("T_ID"))

pM3 <- ggplot(Cluster_Compare_M3_melt, aes(x = T_ID, y = variable)) +geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("M1"="#db5243","M2"="#4db1c7","M3"="#309d88","S-Mb"="#9dc591","S-Pf"="lightgrey","S-Me"="lightgrey",
                               "SI"="#f9c794","SII"="lightgrey","SIII"="lightgrey","S1"="#9492c2","S2"="lightgrey","S3"="lightgrey")) +
  theme_minimal() +labs(x = "Samples", y = "Categories", fill = "Values") +theme(axis.text.x = element_text(angle = 45, hjust = 1))
pM1 + pM2 + pM3  + patchwork::plot_layout(ncol = 1)

#### Oncoplot ####
Select_gene <- c("TP53","CTNNB1","AXIN1","ALB","APOB","RB1","KEAP1","ARID1A","TSC2","KMT2C");
colors <- c("frameshift deletion"="#2E78B3","frameshift insertion"="#A4D38F","nonframeshift deletion" = "#4A9D47","nonframeshift substitution"= "#F59494","stopgain" = "#979FD8","nonsynonymous SNV" = "#4db1c7","splicing" = "#A75C37","NoMut"= "lightgrey")
mutation_list <- lapply(WES_index_sp, function(x) {mutation_matrix <- matrix(0,nrow = length(levels(x$Gene)),ncol = length(levels(x$T_ID)));dim1 <- as.numeric(x$Gene);dim2 <- as.numeric(x$T_ID);for (i in seq_along(dim1)) {mutation_matrix[dim1[i], dim2[i]] <- 1}
rownames(mutation_matrix) <- levels(x$Gene);colnames(mutation_matrix) <- levels(x$T_ID);freq <- rowMeans(mutation_matrix);mutation_matrix <- mutation_matrix[order(freq, decreasing = TRUE), ];mutation_matrix[Select_gene, , drop = FALSE]})%>%.[names(colors)]
alter_fun <- lapply(colors, function(col){function(x,y,w,h) {grid.rect(x,y,w - unit(0.5,"mm"), h*0.8,gp=gpar(fill=col,col=NA))}});alter_fun$background<-function(x, y, w, h) {grid.rect(x, y, w - unit(0.5, "mm"), h * 0.8,gp = gpar(fill = "#e8e7e3", col = NA))}
WES_group <- Clinic_109[,c("T_ID","Cluster")]
mutation_list_M1 <- lapply(mutation_list, function(x){x[,as.character(Clinic_109[,c("T_ID","Cluster")]$T_ID[Clinic_109[,c("T_ID","Cluster")][,2] == "M1"])]})
top_anno <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),Cluster = Clinic_109$Cluster[Clinic_109$Cluster == "M1"],col = list(Cluster = c("M1" = "#db5243")),gp = gpar(col = "white"))
WES_M1 <- oncoPrint(mutation_list_M1,alter_fun = alter_fun,col = colors,row_order = Select_gene,show_column_names = TRUE,top_annotation = top_anno)
WES_M1_Order <-as.character(Clinic_109[,c("T_ID","Cluster")]$T_ID[Clinic_109[,c("T_ID","Cluster")][,2]=="M1"])[WES_M1@column_order]
mutation_list_M2 <- lapply(mutation_list, function(x){x[,as.character(Clinic_109[,c("T_ID","Cluster")]$T_ID[Clinic_109[,c("T_ID","Cluster")][,2] == "M2"])]})
top_anno <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),Cluster = Clinic_109$Cluster[Clinic_109$Cluster == "M2"],col = list(Cluster = c("M2" = "#4db1c7")),gp = gpar(col = "white"))
WES_M2 <- oncoPrint(mutation_list_M2,alter_fun = alter_fun,col = colors,row_order = Select_gene,show_column_names = TRUE,top_annotation = top_anno)
WES_M2_Order <-as.character(Clinic_109[,c("T_ID","Cluster")]$T_ID[Clinic_109[,c("T_ID","Cluster")][,2]=="M2"])[WES_M2@column_order]
mutation_list_M3 <- lapply(mutation_list, function(x){x[,as.character(Clinic_109[,c("T_ID","Cluster")]$T_ID[Clinic_109[,c("T_ID","Cluster")][,2] == "M3"])]})
top_anno <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),Cluster = Clinic_109$Cluster[Clinic_109$Cluster == "M3"],col = list(Cluster = c("M3" = "#309d88")),gp = gpar(col = "white"))
WES_M3 <- oncoPrint(mutation_list_M3,alter_fun = alter_fun,col = colors,row_order = Select_gene,show_column_names = TRUE,top_annotation = top_anno)
WES_M3_Order <-as.character(Clinic_109[,c("T_ID","Cluster")]$T_ID[Clinic_109[,c("T_ID","Cluster")][,2]=="M3"])[WES_M3@column_order]
g1 <- grid.grabExpr(draw(WES_M1, newpage = FALSE));g2 <- grid.grabExpr(draw(WES_M2, newpage = FALSE));g3 <- grid.grabExpr(draw(WES_M3, newpage = FALSE))
WES_Order <- Clinic_109$ID[match(c(WES_M1_Order,WES_M2_Order,WES_M3_Order),Clinic_109$T_ID)]
WES_T_ID <- Clinic_109$T_ID[match(c(WES_M1_Order,WES_M2_Order,WES_M3_Order),Clinic_109$T_ID)]
gridExtra::grid.arrange(g1, g2, g3, ncol = 3)

#### Pathway Heatmap ####
Heatmap_Data <- Data_GSVA %>%set_rownames(.$ID) %>%.[Clinic_109$ID, Path_TvsN_All$Symbol[Path_TvsN_All$FDR < 0.05]] %>%t() %>%as.data.frame() %>%tibble::rownames_to_column("Pathway") %>%mutate(across(-Pathway, as.numeric))
sample.info <- Clinic_109 %>%arrange(Cluster);sample.info[, Heatmap_Data$Pathway] <- t(Heatmap_Data[, sample.info$ID])
Heatmap_Data$P_kruskal.test <- sapply(Heatmap_Data$Pathway, function(sub){kruskal.test(sample.info[[sub]] ~ sample.info$Cluster)$p.value %>% signif(2)})
Heatmap_Data_psig <- Heatmap_Data %>%filter(P_kruskal.test < 0.05) %>%mutate(Catogorize = KEGG_Catogorize$Catogorize[match(Pathway, KEGG_Catogorize$Pathway)]) %>%arrange(Catogorize, P_kruskal.test) %>%tibble::column_to_rownames("Pathway")
y <- as.matrix(Heatmap_Data_psig[, WES_Order]);scale_rows <- function(x){(x - rowMeans(x, na.rm = TRUE)) / apply(x, 1, sd, na.rm = TRUE)};y <- scale_rows(y)
## Annotation
Top_anno = HeatmapAnnotation(
  gp = gpar(col = "white"),
  foo = anno_block(gp = gpar(fill = c("#db5243","#4db1c7","#309d88")),
                   labels = c("M1", "M2", "M3"),labels_gp = gpar(fontsize = 15)),
  na_col = "grey", border = TRUE,show_annotation_name = T)
Right_annotation <- rowAnnotation(pvalue = anno_numeric(Heatmap_Data_psig$P_kruskal.test, rg = c(min(Heatmap_Data_psig$P_kruskal.test), 1),x_convert = function(x) -log10(Heatmap_Data_psig$P_kruskal.test),
                                                        labels_format = function(x) sprintf("%.2e", Heatmap_Data_psig$P_kruskal.test)),annotation_name_rot = 0)
Left_annotation <- rowAnnotation(foo = anno_block(gp=gpar(fill=colors), labels=unique(Heatmap_Data_psig$Catogorize), labels_gp=gpar(fontsize=5)),na_col="grey", border=F, show_annotation_name=F)
####   Draw   ####
ComplexHeatmap::Heatmap(y,col = colorRamp2(c(-2,0,2), c("#377EB8","white","#FF3333")),cluster_rows = FALSE,cluster_columns = FALSE,row_split = Heatmap_Data_psig$Catogorize,column_split = Clinic_109$Cluster[match(c(WES_M1_Order,WES_M2_Order,WES_M3_Order), Clinic_109$T_ID)],
                        top_annotation = Top_anno,left_annotation = Left_annotation,right_annotation = Right_annotation,row_names_gp = gpar(fontsize=15),show_column_names = F,show_row_names = TRUE,column_title = NULL,row_title = NULL,border = NA,
                        heatmap_legend_param = list(title=NULL, legend_height=unit(2,"cm"), legend_direction="vertical"))

####################################
####  Fig1E  #######################
####################################
#### Tyrosine  #######
Tyrosine_Path <-  Meta_109_Raw[,c(6,12:ncol(Meta_109_Raw))]%>%set_rownames(.$ID)%>%.[-1]%>%t()%>%as.data.frame()%>%.[c("C05582","C00042","C01161","C04043","C05587","C03758","C00355","C00082","C00547","C00788","C05588")]%>%
  set_colnames(c("HVA","Succinate","DOPAC","DOPAL","3MT","Dopamine","LDOPA","Tyrosine","Noradrenaline","Adrenaline","Metanephrine"))%>%cbind(sample.name=rownames(.),Cluster=Clinic_109$Cluster[match(rownames(.),Clinic_109$sample.name)],.)%>%.[grepl("", .$Cluster),-1]
means_by_cluster <- Tyrosine_Path %>%group_by(Cluster) %>%summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%mutate(across(-Cluster, ~ scale(.x)[, 1])) %>%arrange(Cluster) %>%mutate(SampleID = paste0("Sample", row_number()))
long_df <- melt(means_by_cluster, id.vars = c("SampleID", "Cluster"));long_df$label <- sprintf("%.2f", long_df$value)
####   Draw   ####
ggplot(long_df, aes(x = SampleID, y = variable, fill = value)) +geom_tile(color = "white") +geom_text(aes(label = label), size = 3) +  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = median(long_df$value, na.rm = TRUE)) +
  facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +labs(title = "All Features Heatmap with Values",x = "Samples", y = "Features") +
  theme_minimal() +theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "grey90", color = NA),strip.text = element_text(size = 10))

####################################
####  Fig1F  #######################
####################################
#### Calculate Cluster Meta ####
Meta_Cluster <- Meta_109_Raw[,c(6,11:ncol(Meta_109_Raw))]%>% set_rownames(.$ID)%>% .[,Sample.info$sample.name[Sample.info$group=="Tumor"]]%>% t()%>% as.data.frame()%>% 
  mutate(sample.name  = rownames(.),Cluster = Clinic_109$Cluster[match(rownames(.),Clinic_109$sample.name)],Type=ifelse(Cluster=="M1","M1","Others"))
Index <- head(colnames(Meta_Cluster), -3)
Meta_Cluster_All <- purrr::map_dfr(Index, function(sub) { Data <- Meta_Cluster %>%dplyr::select(sample.name, Type, FPKM = all_of(sub)) %>%
  filter(!is.na(FPKM)) %>%dplyr::arrange(Type);Other_FPKM <- Data$FPKM[Data$Type == "Others"];M1_FPKM <- Data$FPKM[Data$Type == "M1"]
  tibble(Symbol = sub,Log2FC = log2(mean(M1_FPKM, na.rm = TRUE) / mean(Other_FPKM, na.rm = TRUE)),Pvalue = wilcox.test(Other_FPKM, M1_FPKM)$p.value)})%>%
  mutate(FDR = p.adjust(.$Pvalue,method="fdr"),name = Meta_109_Raw$compound.name[match(Symbol,Meta_109_Raw$ID)],weight_Cluster=sqrt(abs(Log2FC)*(-log2(FDR))))
#### Meta_Clinic ####
Meta_Cli <- Meta_109_Raw[,c(6,11:ncol(Meta_109_Raw))]%>% set_rownames(.$ID)%>% .[,Sample.info$sample.name[Sample.info$group=="Tumor"]]%>% t()%>% as.data.frame()%>% 
  mutate(sample.name  = rownames(.),Overall_survial = Clinic_109$Overall_survial[match(rownames(.),Clinic_109$sample.name)],Survial_status = Clinic_109$Survial_status[match(rownames(.),Clinic_109$sample.name)])
Index1 <- head(colnames(Meta_Cluster), -3)%>%setdiff(.,c("C05835","C01252","C00785"))
Index2 <- c("C05835","C01252","C00785")
Meta_Cli_Risk <- purrr::map_dfr(Index1, function(sub) { Data <- Meta_Cli %>%dplyr::select(sample.name, Survial_status,Overall_survial, FPKM = all_of(sub)) %>%
  filter(!is.na(FPKM))%>%mutate(logFPKM=log2(FPKM +1));Cutpoint=surv_cutpoint(Data,time = "Overall_survial", event="Survial_status", minprop=0.06,variables = "logFPKM");Res.cat=surv_categorize(Cutpoint);
  Res.cat$logFPKM <- factor(Res.cat$logFPKM, levels = c("low", "high"));fit=survfit(Surv(Overall_survial, Survial_status) ~logFPKM, data = Res.cat);Model=survdiff(survival::Surv( Overall_survial, Survial_status) ~  logFPKM,data=Res.cat, na.action=na.exclude);
  KMP= 1-pchisq(Model$chisq, df=length(levels(factor(Res.cat$logFPKM)))-1);cox_model <- coxph(Surv(Overall_survial, Survial_status) ~ logFPKM, data = Res.cat);cox_summary <- summary(cox_model);
  tibble(Symbol = sub,KMP = KMP,HR=cox_summary$coefficients[1, "exp(coef)"],HR_lower=cox_summary$conf.int[1, "lower .95"],HR_upper=cox_summary$conf.int[1, "upper .95"],HR_p=cox_summary$coefficients[1, "Pr(>|z|)"],weight_Cli=sqrt(abs(HR)*(-log2(HR_p))))})
Meta_Cli_Pro <- purrr::map_dfr(Index2, function(sub) { Data <- Meta_Cli %>%dplyr::select(sample.name, Survial_status,Overall_survial, FPKM = all_of(sub)) %>%
  filter(!is.na(FPKM))%>%mutate(logFPKM=log2(FPKM +1));Cutpoint=surv_cutpoint(Data,time = "Overall_survial", event="Survial_status", minprop=0.06,variables = "logFPKM");Res.cat=surv_categorize(Cutpoint);
  Res.cat$logFPKM <- factor(Res.cat$logFPKM, levels = c("high","low"));fit=survfit(Surv(Overall_survial, Survial_status) ~logFPKM, data = Res.cat);Model=survdiff(survival::Surv( Overall_survial, Survial_status) ~  logFPKM,data=Res.cat, na.action=na.exclude);
  KMP= 1-pchisq(Model$chisq, df=length(levels(factor(Res.cat$logFPKM)))-1);cox_model <- coxph(Surv(Overall_survial, Survial_status) ~ logFPKM, data = Res.cat);cox_summary <- summary(cox_model);
  tibble(Symbol = sub,KMP = KMP,HR=cox_summary$coefficients[1, "exp(coef)"],HR_lower=cox_summary$conf.int[1, "lower .95"],HR_upper=cox_summary$conf.int[1, "upper .95"],HR_p=cox_summary$coefficients[1, "Pr(>|z|)"],weight_Cli=sqrt(abs(HR)*(-log2(HR_p))))})
Meta_Cli_All <- rbind(Meta_Cli_Risk,Meta_Cli_Pro);
####   Draw   ####
UpM1 <- filter(Meta_Cluster_All,Log2FC>0&FDR<0.05);Acid <-filter(Meta_Class,Category == "Amino acid");Surv_sig <- filter(Meta_Cli_All,KMP<0.05)
All <- list( UpM1, Acid, Surv_sig) %>%purrr::reduce(inner_join, by = "Symbol")%>%mutate(weight=weight_Cluster+weight_Cli)%>%arrange(desc(weight));setdiff(c("C05835","C01252","C00785"),All[1:10,]$Symbol)
#### venn ####
venn_list <- list(`Upregulated in M1 (P < 0.05)` = UpM1$Symbol,`Prognostic in cohort` = Surv_sig$Symbol,`Amino acids` = Acid$Symbol)
ggvenn::ggvenn(venn_list,fill_color = c("#4db1c7","#db5243","#309d88"),stroke_size = 0.7,set_name_size = 4,text_size = 4,fill_alpha = 0.5)
#### Weight ####
ggplot(All[1:10,], aes(x = reorder(name.x, weight), y = weight)) +geom_segment(aes(x = name.x, xend = name.x, y = 1, yend = weight), color = "gray") +geom_point(aes(size = weight), color = "#ef8c62") +
  scale_size_continuous(range = c(2, 6)) +coord_flip() +scale_y_continuous(limits = c(1, NA)) +theme_minimal() +labs(x = NULL, y = "Frequency", title = "Lollipop Chart (Size by Frequency)") +theme(axis.text.y = element_text(size = 10))



