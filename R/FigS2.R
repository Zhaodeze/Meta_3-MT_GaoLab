
rm(list = ls())

Target_3MT <-  data.table::fread("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/All_Files/Mice_3MT.csv")%>% 
  as.data.frame()

# identical(sample.Tumor$Sample_ID, sample.Tumor$ID)  # 严格比较，类型和值都一致才返回 TRUE
CD8 <- data.table::fread("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/All_Files/IHC/mIHC_Mice_20260108.csv",header = T)  %>% 
  as.data.frame()%>% 
  cbind(Order=c(1:nrow(.)),Target_3MT=Target_3MT$nmol_g[match(.$Sample_ID,Target_3MT$Sample_ID)])
cols_to_convert <- sapply(CD8, function(x) all(grepl("%$", x)))
CD8[cols_to_convert] <- lapply(CD8[cols_to_convert], function(x) as.numeric(sub("%", "", x)))# / 100)
sem <- function(x, na.rm = TRUE) {
  x <- if (na.rm) x[!is.na(x)] else x
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}
# colnames(CD8_stats)

CD8_stats <- CD8 %>%
  group_by(Sample_id, Sample_ID, Target_3MT) %>%
  summarise(
    n = sum(!is.na(Single_CD8)),                    # 或 n()，看你是否有NA
    Single_CD8_m   = mean(Single_CD8, na.rm = TRUE),
    Single_CD8_sem = sem(Single_CD8),
    Double_Positive_m   = mean(Double_Positive, na.rm = TRUE),
    Double_Positive_sem = sem(Double_Positive),
    .groups = "drop"
  ) %>%
  arrange(desc(Target_3MT), desc(Single_CD8_m)) %>%
  .[!.$Sample_id %in% c("1MYC_sgP53_1","1MYC_sgP53_3",
                        "2MYC_PI3K_1","2MYC_PI3K_4","2MYC_PI3K_7",
                        "3CON MC-2",
                        "4CTRL3",
                        "5MYC_YAP_2",
                        "6PI3K_CCND1_2","6PI3K_CCND1_6",
                        "8CTNNB1_MYC_1","8CTNNB1_MYC_6",
                        "9AKT_MET_4",
                        "10HEPA1-6 1",
                        "12CTNNB1_AKT_1","12CTNNB1_AKT_6",
                        "13AKT_NRAS_2","13AKT_NRAS_4"
  ),] %>%
  cbind(Group = ifelse(.$Target_3MT >= median(.$Target_3MT), "High", "Low"), .)

CD8_stats_sem <- CD8_stats %>%
  group_by(Sample_ID, Target_3MT, Group) %>%
  summarise(
    n_rep = n(),
    Single_CD8_mean = mean(Single_CD8_m, na.rm = TRUE),
    Single_CD8_sem  = sd(Single_CD8_m, na.rm = TRUE) / sqrt(n_rep),
    Double_Positive_mean = mean(Double_Positive_m, na.rm = TRUE),
    Double_Positive_sem  = sd(Double_Positive_m, na.rm = TRUE) / sqrt(n_rep),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_sem"), ~replace_na(., 0))) %>%
  arrange(desc(Target_3MT))


CD8_stats_sd$Sample_ID <- factor(CD8_stats_sd$Sample_ID ,levels = CD8_stats_sd$Sample_ID )
sample_order<- CD8_stats_sd$Sample_ID    

colors<-  c(
  "#DC050C", "#FB8072", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7")


theme_clean_xy <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    
    axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p1 <- ggplot(CD8_stats_sem,aes(x = Target_3MT, y = Single_CD8_mean)) +
      geom_point(aes(color = Sample_ID, size = Single_CD8_mean),alpha = 0.9) +
      scale_size_continuous(name = "Mean\n(Single CD8)",range = c(2.5, 7)) +
      ggnewscale::new_scale("size") +
      geom_point(
      aes(color = Sample_ID, size = Single_CD8_sem),shape = 1,stroke = 0.6) +
      scale_size_continuous(name = "SEM\n(Single CD8)",range = c(2, 10)) +
      geom_smooth(method = "lm", color = "#594a41", se = TRUE) +
      stat_cor(method = "pearson") +
      scale_color_manual(values = setNames(rev(colors), sample_order),name = "Sample ID") +
      labs(x = "Meta_3MT",y = "Single CD8",
           title = "Correlation between Meta_3MT and Single CD8",
           subtitle = "Solid dot = mean, open circle = SEM") +
      theme_clean_xy +
       guides(color = guide_legend(order = 1),size  = guide_legend(order = 2))



p2 <- ggplot(CD8_stats_sem,aes(x = Target_3MT, y = Double_Positive_mean)) +
  geom_point(aes(color = Sample_ID, size = Double_Positive_mean),alpha = 0.9) +
  scale_size_continuous(name = "Mean\n(GZMB CD8)",range = c(2.5, 7)) +
  ggnewscale::new_scale("size") +
  geom_point(
    aes(color = Sample_ID, size = Single_CD8_sem),shape = 1,stroke = 0.6) +
  scale_size_continuous(name = "SEM\n(GZMB CD8)",range = c(2, 10)) +
  geom_smooth(method = "lm", color = "#594a41", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = setNames(rev(colors), sample_order),name = "Sample ID") +
  labs(x = "Meta_3MT",y = "GZMB CD8",
       title = "Correlation between Meta_3MT and GZMB CD8",
       subtitle = "Solid dot = mean, open circle = SEM") +
  theme_clean_xy +
  guides(color = guide_legend(order = 1),size  = guide_legend(order = 2))



p1 + p2   + patchwork::plot_layout(ncol = 2)

ggsave("/home/data/user/shenxia/HCC/met_subtype/0_New_results/MIce/mIHC_Mice_Cor_CD8_Targets.pdf", width = 16, height = 6)








rm(list = ls())
colors<-  c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig1/Fig2_All_DataMatrix.Rdata")

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


Mice_3MT <-  data.table::fread("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/All_Files/Mice_3MT.csv")

Mice_3MT$Name = factor(Mice_3MT$Sample_ID,levels = rev(Mice_3MT$Sample_ID))
Mice_3MT$Value <- round(Mice_3MT$nmol_g,2)

ggplot(Mice_3MT, aes(x = Value, y = Name)) +
  geom_point(aes(color = Name), size = 8) +
  geom_point(aes(color = Name), size = 10, shape = 21, fill = NA) +
  geom_segment(aes(x = 0, xend = Value, y = Name, yend = Name, color = Name), linewidth = 1) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.position = "none",  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5), 
        panel.grid.minor = element_line(color = "grey", linewidth = 0.25)) + 
  xlim(0, 30)+
  # geom_text(aes(label = Value), hjust = -1) +
  labs(x = "3-MT(nmol/g)", y = "") 

# ggsave("/home/data/user/shenxia/HCC/met_subtype/0_New_results/MIce/lollipop_Target_3MT.pdf", width = 6, height = 6)
