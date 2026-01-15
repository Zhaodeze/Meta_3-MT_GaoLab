library(ggplot2);library(dplyr);library(magrittr);library(tidyr);
library(changepoint);library(ggplot2);library(ggh4x);library(patchwork);

rm(list = ls())
load("/home/data/user/shenxia/HCC/met_subtype/0_New_results/Fig7/Fig7_All_DataMatrix.Rdata")
source("/home/data/user/shenxia/HCC/met_subtype/00_Code_deze/0_Data/Function/00_Functions_used.R")

####################################
####  Fig7 B  ######################
####################################

ggplot() +
  geom_segment(data = AllTimes,aes(x = Start, xend = End,
               y = Patient, yend = Patient, color = Outcome),linewidth = 4) +
  scale_color_manual(values = c("PD" = "#cbbea7","PR" = "#88bbbd","SD" = "#bcbdc0")) +
  geom_point(data = Points,aes(x = Time, y = Patient, shape = Outcome),size = 3) +
  scale_shape_manual(values =  c("PD" = 16,"PR" = 17,"SD" = 15 )) +
  labs(x = "Time (days)", y = "Patient",shape = "Treatment Outcome") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor = element_blank())+
  theme_classic()

####################################
####  Fig7 C  ######################
####################################
ggplot(T_Change_RECIST1.1, aes(x = Patients, y = Change, fill = Outcome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("PD" = "#cbbea7","PR" = "#88bbbd","SD" = "#bcbdc0")) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#cbbea7")+
  geom_hline(yintercept = -30, linetype = "dashed", color = "#88bbbd")+
  theme_bw() +
  labs(y = "Change from baseline (%)", x = "Patient ID",title = "") +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank())

ggplot(T_Change_mRECIST, aes(x = Patients, y = Change, fill = Outcome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("PD" = "#cbbea7","PR" = "#88bbbd","SD" = "#bcbdc0")) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#cbbea7")+
  geom_hline(yintercept = -30, linetype = "dashed", color = "#88bbbd")+
  theme_bw() +
  labs(y = "Change from baseline (%)", x = "Patient ID",title = "") +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank())
####################################
####  Fig7 GL ######################
####################################
p1 <- ggplot(Clinic_Tumor_Change,
             aes(x = Time, y = Value,
                 group = PatientID,
                 color = PatientID)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Group, scales = "free_y",axes = "all") +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Tumor trajectories") +
  theme_clean_xy +
  theme(legend.position = "none")

p2 <- ggplot(Clinic_3MT_Change,
             aes(x = Time, y = Mean, group = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_line(linewidth = 1, color = "grey40") +
  geom_point(size = 2, color = "red") +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SEM),
        ymax = Mean + SEM),
    width = 0.15,
    color = "pink"
  ) +
  facet_wrap(~ Group, scales = "free_y",axes = "all") +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(title = "3MT") +
  theme_clean_xy

p3 <- ggplot(Clinic_CD8_Change,
             aes(x = Time, y = Value, group = PatientID)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_line(linewidth = 1, color = "grey") +
  geom_point(size = 2, color = "darkblue") +
  facet_wrap(~ PatientID, scales = "free_y",axes = "all") +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "CD8 trajectories") +
  theme_clean_xy
final_plot <- p1 | p2 | p3
