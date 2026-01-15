library(tidyr);library(dplyr);library(ggplot2);library(ggpubr);library(magrittr);library(circlize);library(purrr);library(ggrepel)


rm(list = ls())
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig3/Fig3_All_DataMatrix.Rdata")
source("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/Function/00_Functions_used.R")
####################################
####  FigS3F  ######################
####################################
NoSig_points<-mRNA_Data%>%filter(Mark!="yes");Mark_points<-mRNA_Data%>%filter(Mark=="yes");
label_points<-mRNA_Data%>%filter(Mark=="yes"&LogP>-log10(0.05)&abs(log2FoldChange)>1)
Colors <- c("no" ='grey',  "Nosig" ='grey', "yes"='#603813')

ggplot(mRNA_Data, aes(x = log2FoldChange, y = LogP)) +
  geom_point(data = NoSig_points, aes(color = Mark, size = abs(log2FoldChange)), alpha = 0.6) +
  geom_point(data = Mark_points, aes(color = Mark, size = abs(log2FoldChange)), alpha = 0.6) +
  scale_color_manual(values = Colors) +theme_bw(base_size = 16) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "Moccasin", lwd = 1) +labs(x = "LogFC", y = "-Log10P") +
  geom_text_repel(data = label_points,aes(label = Symbol),size = 5,color = "black",segment.color = "black",segment.size = 0.5,
                  segment.curvature = 0.1,segment.angle = 20,segment.ncp = 3,box.padding = 0.6,point.padding = 0.4,min.segment.length = 0)+
  theme_clean_xy

