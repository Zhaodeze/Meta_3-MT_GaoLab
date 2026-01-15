
rm(list = ls())
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig3/Fig3_All_DataMatrix.Rdata")
source("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/Function/00_Functions_used.R")
####################################
####  Fig3 A  ######################
####################################

UP_GO_summary <- UP_GO %>%group_by(Class) %>%summarise(mean_score = mean(abs(Enrichment), na.rm = TRUE),n_pathways = n()) %>%arrange(desc(mean_score))
ggplot(UP_GO_summary, aes(x = reorder(Class, mean_score), y = mean_score)) +geom_col(fill = "steelblue") +coord_flip() +
  labs(y = "Mean |Enrichment Score|", x = "Class",title = "Pathway Category Enrichment Trend Up")

Down_GO_summary <- Down_GO %>%group_by(Class) %>%summarise(mean_score = mean(abs(Enrichment), na.rm = TRUE),n_pathways = n()) %>%arrange(desc(mean_score))
ggplot(Down_GO_summary, aes(x = reorder(Class, mean_score), y = mean_score)) +geom_col(fill = "steelblue") +coord_flip() +
  labs(y = "Mean |Enrichment Score|", x = "Class",title = "Pathway Category Enrichment Trend Down")

####################################
####  Fig3 B  ######################
####################################
ggplot(Impact_Pathway,aes(x = Impact,y = `LOG10(p)`,color = group)) +geom_point(size = 4, alpha = 0.85) +
  labs(x = "Impact",y = "Log10(p)",title = "Pathway Impact vs Significance") +
  theme_classic(base_size = 18)

####################################
####  Fig3 H  ######################
####################################
lim_expr <- max(Arg_mRNA$MeanExpr, na.rm = TRUE)
ggplot(Arg_mRNA,aes(x = Group, y = Symbol, fill = MeanExpr)) +
  geom_tile(width = 0.92, height = 0.92,color = "white", linewidth = 0.6) +
  geom_text(aes(label = round(MeanExpr,2) , color = txt_col),size = 4, fontface = "bold") +
  scale_color_identity() +
  geom_text(data = subset(Arg_mRNA, Group == "Case"),aes(label = Star),
    position = position_nudge(x = 0.30, y = -0.30),size = 5,fontface = "bold") +
  scale_fill_gradient2(low = "#2b6cb0",mid = "#f7fafc",high = "#c53030",
    midpoint = 0.58,limits = c(0, lim_expr),oob = squish,name = "log2(Case/Ctrl)") +
  coord_fixed() +
  labs(x = NULL, y = NULL,title = "Arginine metabolism") +
  theme_classic(base_size = 13) 
  

