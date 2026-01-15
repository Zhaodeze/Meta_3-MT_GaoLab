####################################
####  Fig2 Prepare #################
####################################
rm(list = ls())
library(tidyr);library(dplyr);library(magrittr);library(data.table);library(xCell);library(tidyverse);library(ggpubr);library(tibble);library(readr);
library(purrr);library(circlize);library(ggplot2);library(ggpubr)  
library(grid)
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig2/Fig2_All_DataMatrix.Rdata")
source("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/Function/00_Functions_used.R")
####################################
####  Fig2A  #######################
####################################

####################################
#### Data Prepare #######
# Result_28 <- BayesDeBulk(n.iter = 10000, burn.in = 1000,
#                          Y = list(Expr_mRNA,Expr_Protein),
#                          markers = Signature28_Markers)
# Result_Bindea <- BayesDeBulk(n.iter = 10000, burn.in = 1000,
#                              Y = list(Expr_mRNA,Expr_Protein),
#                              markers = Signature_Bindea_Markers)
####################################
#### Cor  #######
Immune28_data <- Result_28$cell.fraction%>%as.data.frame()%>%
  cbind(T_ID=rownames(.),Log_3MT=Clinic_109$Meta_3MT[match(rownames(.),Clinic_109$T_ID)]%>%log2(),.)%>%
  mutate(across(where(~ all(!is.na(suppressWarnings(as.numeric(.))))),   ~ as.numeric(.)))%>%
  arrange(Log_3MT)
ImmuneBindea_data <- Result_Bindea$cell.fraction%>%as.data.frame()%>%
  cbind(T_ID=rownames(.),.)%>%
  mutate(across(where(~ all(!is.na(suppressWarnings(as.numeric(.))))),   ~ as.numeric(.)))
TME_109 <- Immune28_data %>%inner_join(ImmuneBindea_data, by= "T_ID")%>%set_rownames(.$T_ID)%>%dplyr::select(-T_ID)

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

####################################
####  Fig2B    #####################
####################################
Order <- arrange(Mice_3MT,Value)
ggplot(Mice_3MT, aes(x = Value, y = Sample_ID)) +
  geom_point(aes(color = Sample_ID), size = 8) +
  geom_point(aes(color = Sample_ID), size = 10, shape = 21, fill = NA) +
  geom_segment(aes(x = 0, xend = Value, y = Sample_ID, yend = Sample_ID, color = Sample_ID), linewidth = 1) +
  xlim(0, 30)+
  labs(x = "3-MT(nmol/g)", y = "") +
  scale_y_discrete(limits = Order$Sample_ID) +
  scale_color_manual(values = colors) + 
  theme_classic() +
  theme(legend.position = "none",  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5), 
        panel.grid.minor = element_line(color = "grey", linewidth = 0.25))

####################################
####  Fig2C    #####################
####################################
p1 <- ggplot(CD8_stats_sem,aes(x = Mice_3MT, y = Single_CD8_mean)) +
  geom_point(aes(color = Sample_ID, size = Single_CD8_mean),alpha = 0.9) +
  scale_size_continuous(name = "Mean\n(Single CD8)",range = c(1, 5)) +
  ggnewscale::new_scale("size") +
  geom_point(aes(color = Sample_ID, size = Single_CD8_mean+Single_CD8_sem),shape = 1,stroke = 0.6) +
  scale_size_continuous(name = "SEM\n(Single CD8)",range = c(1, 7)) +
  geom_smooth(method = "lm", color = "#594a41", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = colors) + 
  labs(x = "Meta_3MT",y = "Single CD8",
       title = "Correlation between Meta_3MT and Single CD8",
       subtitle = "Solid dot = mean, open circle = SEM") +
  theme_clean_xy +
  guides(color = guide_legend(order = 1),size  = guide_legend(order = 2))



p2 <- ggplot(CD8_stats_sem,aes(x = Mice_3MT, y = Double_Positive_mean)) +
  geom_point(aes(color = Sample_ID, size = Double_Positive_mean),alpha = 0.9) +
  scale_size_continuous(name = "Mean\n(GZMB CD8)",range = c(1, 5)) +
  ggnewscale::new_scale("size") +
  geom_point(aes(color = Sample_ID, size = Double_Positive_sem),shape = 1,stroke = 0.6) +
  scale_size_continuous(name = "SEM\n(GZMB CD8)",range = c(1, 7)) +
  geom_smooth(method = "lm", color = "#594a41", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = colors) + 
  labs(x = "Meta_3MT",y = "GZMB CD8",
       title = "Correlation between Meta_3MT and GZMB CD8",
       subtitle = "Solid dot = mean, open circle = SEM") +
  theme_clean_xy +
  guides(color = guide_legend(order = 1),size  = guide_legend(order = 2))



p1 + p2   + patchwork::plot_layout(ncol = 2)



