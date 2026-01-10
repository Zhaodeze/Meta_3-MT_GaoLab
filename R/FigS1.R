####################################
####  Fig1 Prepare #################
####################################
rm(list = ls())
library(ggplot2);library(dplyr);library(magrittr);library(tidyr);library(matrixStats);library(circlize);library(GSVA);#ver.1.40.1
library(RColorBrewer);library(ComplexHeatmap);library(survival);library(reshape2);library(survminer);library(plotly);library(ropls);library(NMF);
#Fig1A  Meta_TvsN_All
#Fig1B  Impact_Pathway
#Fig1C  sample.TN,p_hsa,Meta_109,Clinic_109,Data_GSVA,KEGG_Catogorize
#Fig1D  Cluster_Compare,WES_index_sp,Path_TvsN_All,WES_Order
#Fig1E  Meta_109_Raw
#Fig1F  Sample.info,Meta_Class
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig1/Fig1_All_DataMatrix.Rdata")
####################################
####  FigS1A  #######################
####################################
Pieplot <- table(Meta_TvsN_All$Category)%>%as.data.frame()%>%
  mutate(perc=Freq/sum(Freq)*100,label=paste0(Var1," (",round(perc,3),"%)"),
         Category = factor(Var1,levels = c("Amino acid", "Carbohydrates","Lipid", "Nucleotide","Vitamins","Peptide", "Xenobiotics","Other")))
ggplot(Pieplot, aes(x = "", y = perc, fill = Category)) +geom_bar(stat = "identity", width = 1, color = "white") +coord_polar(theta = "y") +
  theme_void(base_size = 14) +scale_fill_manual(values = c('#81cac0','#b1b0d0','#f37d74','#77a8c6','#f3c2da','#a6d177','#ac7cb7','#f5a86a')) +
  geom_text(aes(label = paste0(round(perc, 1), "%")),position = position_stack(vjust = 0.5),color = "black", size = 4) +labs(fill = "Category")



####################################
####  FigS1B-D  ####################
####################################
data_pca <- select(Meta_109_Raw,Sample.info$sample.name)%>% apply(.,2,function(x)log2(x+1))
final_data <- as.data.frame(prcomp(t(data_pca), scale. = T)$x)%>%
              mutate(Group=Sample.info$group[match(rownames(.),Sample.info$sample.name)])
fig <- plot_ly(final_data, x = ~PC1, y = ~PC2, z = ~PC3,
               color = ~final_data$Group, colors = c('#0e7176','#e69f00','#fe7801'),
               marker = list(size = 6,line = list(width = 1,color = 'rgb(230,230,230)'))) %>%add_markers(size = 22)
Explained_variance_ratio <- 100 * sum(summary(prcomp(t(data_pca), scale. = T))[["importance"]]['Proportion of Variance',])
fig <- fig %>%layout(title = paste0('Total Explained Variance = ',Explained_variance_ratio),scene = list(bgcolor = "white"))
fig



data_pca <- select(Meta_109_Raw,Sample.info$sample.name[Sample.info$group!="QC"])%>% apply(.,2,function(x)log2(x+1))
x <- t(data_pca);y <- factor(Sample.info$group[match(rownames(x),Sample.info$sample.name)]);
oplsda_model <- opls(x, y, predI = 1, orthoI = NA)
plot(oplsda_model, type = "permutation")
plot_df <- data.frame(t1 = oplsda_model@scoreMN[,1],o1 = oplsda_model@orthoScoreMN[,1],o2 = oplsda_model@orthoScoreMN[,2],group = y)
fig <- plot_ly(plot_df, x = ~t1, y = ~o1, z = ~o2,color = ~plot_df$group, colors = c('#e69f00','#0e7176'),
               marker = list(size = 6,line = list(width = 1,color = 'rgb(230,230,230)'))) %>%add_markers(size = 22)
tit = c(paste0("R2Y=",oplsda_model@summaryDF$`R2Y(cum)`," ","Q2=",oplsda_model@summaryDF$`Q2(cum)`))
fig <- fig %>%layout(title = tit,scene = list(bgcolor = "white"))
fig

####################################
####  FigS1E  ####################
####################################

Path_TvsN_All$Catogorize <- factor(Path_TvsN_All$Catogorize, levels =  c("Amino acid","Drug","Carbohydrate","Lipid","Other"))
Path_TvsN_All <- arrange(Path_TvsN_All,desc(Catogorize),Difference)
Path_TvsN_All$Symbol <- factor(Path_TvsN_All$Symbol, levels = Path_TvsN_All$Symbol)
Path_TvsN_All$log10FDR <- -log10(Path_TvsN_All$FDR)

Path_TvsN_Sig <- filter(Path_TvsN_All,FDR<0.05)
Path_TvsN_Sig$Catogorize <- paste0(Path_TvsN_Sig$Catogorize,"1")
Point_color_1 <- c("Other1"  ='#e4a91f',"Lipid1"='#2d8756',"Carbohydrate1"="#0070b8","Drug1"='#c99f62',"Amino acid1"='#ea5517')
Point_color_1_light <- sapply(Point_color_1, function(color) alpha(color, 0.3)) 
names(Point_color_1_light) <- gsub("1","",names(Point_color_1_light))
Point_color <- c(Point_color_1,Point_color_1_light)

ggplot() +  geom_segment(data = Path_TvsN_All,aes(x = 0, xend = Difference, y = Symbol, yend = Symbol,color = Catogorize), linewidth = 0.5) +
  geom_point(data = Path_TvsN_Sig,aes(x = Difference, y = Symbol, color = Catogorize, size = log10FDR+2), shape = 1) + 
  geom_point(data = Path_TvsN_All,aes(x = Difference, y = Symbol,color = Catogorize, size = log10FDR), shape = 16) +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +  
  scale_color_manual(values = Point_color) +theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
        axis.line = element_line(color = "black"), 
        panel.grid.major = element_line(color = "grey", linewidth = 0.25), 
        panel.grid.minor = element_line(color = "grey", linewidth = 0.25)) 

####################################
####  Fig1F  #######################
####################################

#### Arginine  #######
Arginine_Meta <- intersect(p_hsa$`Arginine and proline metabolism - Homo sapiens (human)`,Meta_109_Raw$ID)
Arginine_Path <-  Meta_109_Raw[,c(6,10,12:ncol(Meta_109_Raw))]%>%filter(ID%in%Arginine_Meta)%>%set_rownames(.$compound.name)%>%.[-c(1,2)]%>%t()%>%as.data.frame()%>%cbind(sample.name=rownames(.),Cluster=Clinic_109$Cluster[match(rownames(.),Clinic_109$sample.name)],.)%>%
  .[grepl("", .$Cluster),-1]
means_by_cluster <- Arginine_Path %>%group_by(Cluster) %>%summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%mutate(across(-Cluster, ~ scale(.x)[, 1])) %>%arrange(Cluster) %>%mutate(SampleID = paste0("Sample", row_number()))
long_df <- melt(means_by_cluster, id.vars = c("SampleID", "Cluster"));long_df$label <- sprintf("%.2f", long_df$value)
####   Draw   ####
ggplot(long_df, aes(x = SampleID, y = variable, fill = value)) +geom_tile(color = "white") +geom_text(aes(label = label), size = 3) +  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = median(long_df$value, na.rm = TRUE)) +
  facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +labs(title = "All Features Heatmap with Values",x = "Samples", y = "Features") +
  theme_minimal() +theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "grey90", color = NA),strip.text = element_text(size = 10))


