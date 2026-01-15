library(tidyr);library(dplyr);library(ggplot2);library(ggpubr);library(magrittr);
library(circlize);library(purrr);library(ggrepel)

rm(list = ls())
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig4/Fig4_All_DataMatrix.Rdata")
source("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/Function/00_Functions_used.R")

####################################
####  Fig4 B  ######################
####################################

ggplot(Mass_Data, aes(x = LogUnique_Peptides, y = LogAbundances)) +
  geom_point(aes(color = Change,size=Unique_Peptides), alpha = 0.6) +
  scale_color_manual(values = Colors) +
  geom_vline(xintercept = log2(1.8),lty = 4, col = "Moccasin", lwd = 1) +
  geom_hline(yintercept = log2(2),lty = 4, col = "Moccasin", lwd = 1) + #
  labs(x = "Log2(Unique Peptides)",y = "Log2(Abundances Ratio)") +
  xlim(0, 6) +ylim(-4,7)+
  geom_segment(data = label_points,aes(x = LogUnique_Peptides, xend = LogUnique_Peptides-0.1, y = LogAbundances, yend = LogAbundances+0.1),
               color = "black", linetype = "dashed") +
  geom_text(data = label_points,aes(x = LogUnique_Peptides-0.1, y = LogAbundances+0.1, label = Description),
            size = 5, hjust = 0, vjust = 0, color = "black")+
  theme_classic(base_size = 13) +
  theme_clean_xy


