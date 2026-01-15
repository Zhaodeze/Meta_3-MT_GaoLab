library(ggplot2);library(dplyr);library(magrittr);library(tidyr);
library(changepoint);library(ggplot2);library(ggh4x);library(patchwork);

rm(list = ls())

load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig7/Fig7_All_DataMatrix.Rdata")

####################################
####  FigS7 A ######################
####################################
HCC_data <- Plasma_3MT %>% filter(Group %in% c("HCC"))
N_data <- Plasma_3MT %>% filter(Group %in% c("Normal"))
cp <- cpt.meanvar(sort(N_data$Meta_3MT),method = "PELT",penalty = "MBIC")
Changepoint  <-  tail(sort(N_data$Meta_3MT)[cp@cpts], 2)[1]
df_plot <- rbind(
  data.frame(value = N_data$Meta_3MT_scale, Group = "Normal"),
  data.frame(value = HCC_data$Meta_3MT_scale, Group = "HCC"))
ggplot(df_plot, aes(x = value, fill = Group)) +geom_density(alpha = 0.35) +
  geom_vline(xintercept = log10(Changepoint), linetype = "dashed", color = "red")+
  annotate("text", x = log10(Changepoint)+0.1, y = 6, label =paste0("Threshold=",Changepoint) , hjust = 0, size = 4, color = "red") +
  geom_jitter(data = filter(df_plot, Group == "Normal"),aes(x = value, y =3),
              width = 0, height = 2, alpha = 0.6, color = "grey70", size = 2) +
  scale_fill_manual(values = c("Normal" = "grey70","HCC" = "#d95f02")) +
  scale_x_continuous(breaks = log10(c(0.5, 1, 2, 5, 10)),labels = c("0.5", "1", "2", "5", "10"))+
  labs(x = "Plasma 3-MT",y = "Density",title = "Reference distribution of plasma 3-MT") +
  theme_classic()
