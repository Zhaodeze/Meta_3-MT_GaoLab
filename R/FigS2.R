
rm(list = ls())
####################################
####  FigS2E  ######################
####################################
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig2/Fig2_All_DataMatrix.Rdata")

neg_sig_count <- Immune_Meta_Cor_all %>%
  filter(sig, direction == "neg") %>%
  dplyr::count(Gene, name = "NegSig_n")

df_sig <- Immune_Meta_Cor_all %>% filter(p < 0.05) 
df_nonsig <- Immune_Meta_Cor_all %>% filter(p >= 0.05) 

# 绘图
p <-  ggplot() +
  geom_point(data = df_nonsig, aes(x = Gene, y = Cancer, size = logp, fill = cor),
             shape = 21, color = "white", stroke = 0.3) +
  geom_point(data = df_sig, aes(x = Gene, y = Cancer, size = logp, fill = cor),
             shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Corr") +
  scale_size(range = c(3, 8), name = expression(-log[10](p))) +
  scale_x_discrete(limits = gene_order)  +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )+
  labs(title = "Cor Between 3-MT and Immune28")

####################################
####  FigS2F  ######################
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