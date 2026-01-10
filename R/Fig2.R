####################################
####  Fig2 Prepare #################
####################################
rm(list = ls())
library(tidyr);library(dplyr);library(magrittr);library(data.table);library(xCell);library(tidyverse);library(ggpubr);library(tibble);library(readr);
library(purrr);library(circlize);library(ggplot2);library(ggpubr)  
library(grid)
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig1/Fig2_All_DataMatrix.Rdata")
source("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/Function/00_Functions_used.R")
####################################
####  Fig2A  #######################
####################################

####################################
#### Data Prepare #######
# ImmuneMarkers28 <- Signature_to_markers_no_self(Signature28, names(Signature28))
# Result_28 <- BayesDeBulk(n.iter = 10000, burn.in = 1000,
#                          Y = list(Expr_mRNA,Expr_Protein),
#                          markers = ImmuneMarkers28)
# ImmuneMarkers29 <- Signature_to_markers_no_self(Signature29, names(Signature29))
# Result_29 <- BayesDeBulk(n.iter = 10000, burn.in = 1000,
#                          Y = list(Expr_mRNA,Expr_Protein),
#                          markers = ImmuneMarkers29)

####################################
#### Cor  #######
# Immune28_data <- Result_28$cell.fraction%>%as.data.frame()%>%
#   cbind(T_ID=rownames(.),Log_3MT=Clinic_109$Meta_3MT[match(rownames(.),Clinic_109$T_ID)]%>%log2(),.)%>%
#   mutate(across(where(~ all(!is.na(suppressWarnings(as.numeric(.))))),   ~ as.numeric(.)))%>%
#   arrange(Log_3MT)
# Immune29_data <- Result_29$cell.fraction%>%as.data.frame()%>%
#   cbind(T_ID=rownames(.),.)%>%
#   mutate(across(where(~ all(!is.na(suppressWarnings(as.numeric(.))))),   ~ as.numeric(.)))
# TME_109 <- Immune28_data %>%inner_join(Immune29_data, by= "T_ID")%>%set_rownames(.$T_ID)%>%dplyr::select(-T_ID)

####################################
#### Circos  #######
Cor_Data_Sig <- Cor_Data_Sig%>% arrange(Correlation) 
Heatmap_Data <- TME_109[,c("Log_3MT",as.character(Cor_Data_Sig$Variable))]%>%t()%>%scale_rows(.)%>%t()



circos.clear()
circos.par(start.degree = -30, gap.degree = 120, points.overflow.warning = FALSE)
circos.heatmap(
  Heatmap_Data,
  col = colorRamp2(breaks = c(-1, 0, 1),colors = c("#B3A9EB", "white", "#ffa500")),
  cell.border = "white",
  cell.lwd = 0.1,
  track.height = 0.6,
  cluster = FALSE,
  dend.side = "inside",
  dend.track.height = 0.18,
  rownames.cex = 0.001,
  split = NULL )
ggplot(Cor_Data_Sig, aes(x = Target, y = Variable, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#FF3333", midpoint = 0,
                       limits = c(-max(abs(Cor_Data_Sig$Correlation), na.rm = TRUE), 
                                   max(abs(Cor_Data_Sig$Correlation), na.rm = TRUE))) +
  scale_y_discrete(limits = rev(unique(Cor_Data_Sig$Variable)))+
  geom_text(aes(label = Sig), color = "black", size = 5) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank())
ggsave("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig2/Fig2A_Cor_Cor_Red.pdf", width = 6, height = 6)

####################################
#### Dotplot  #######
Cor_Data <- TME_109%>%as.data.frame()%>%
  .[c("Log_3MT","Activated_CD8_T_cell_28","Effector cell traffic_29")]%>%
  mutate(Cluster = Clinic_109$Cluster[match(rownames(.),Clinic_109$T_ID)],
         M1Other = dplyr::case_when(Cluster == "M1" ~ "M1",
                                    TRUE ~ "Others"))
p_CD8 <-ggplot(Cor_Data, aes(x = Log_3MT, y = `Activated_CD8_T_cell_28`)) +
  geom_point(aes(color = M1Other), alpha = 0.85) +
  stat_ellipse(aes(fill = M1Other),geom = "polygon",type = "norm",
               level = 0.95,alpha = 0.18,color = NA) +
  scale_fill_manual(values = c("M1" = "#db5243", "Others" = "#4db1c7")) +
  scale_color_manual(values = c("M1" = "#db5243", "Others" = "#4db1c7")) +
  geom_smooth(method = "lm", se = TRUE, color = "#594a41", size = 0.8) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top",
           size = 5,color = "black",label.sep = "\n",cor.coef.name = "pearson") +
  theme_minimal(base_size = 14) +
  labs(title = "",x = "Log2(3MT)",y = "Activated CD8 T cell score")
ggExtra::ggMarginal(p_CD8,type = "density",groupColour = TRUE,groupFill = TRUE,alpha = 0.6)

p_Effec <-ggplot(Cor_Data, aes(x = Log_3MT, y = `Effector cell traffic_29`)) +
  geom_point(aes(color = M1Other), alpha = 0.85) +
  stat_ellipse(aes(fill = M1Other),geom = "polygon",type = "norm",
               level = 0.95,alpha = 0.18,color = NA)+
  scale_fill_manual(values = c("M1" = "#db5243", "Others" = "#4db1c7")) +
  scale_color_manual(values = c("M1" = "#db5243", "Others" = "#4db1c7")) +
  geom_smooth(method = "lm", se = TRUE, color = "#594a41", size = 0.8) +
  stat_cor(method = "pearson",label.x.npc = "left",label.y.npc = "top",
           size = 5,color = "black",label.sep = "\n",cor.coef.name = "pearson") +
  theme_minimal(base_size = 14) +
  labs(title = "",x = "Log2(3MT)",y = "Effector cell traffic score")
ggExtra::ggMarginal(p_Effec,type = "density",groupColour = TRUE,groupFill = TRUE,alpha = 0.6)






