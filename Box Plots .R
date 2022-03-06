

library(survival)
library(survminer)
library(dplyr)
library(survminer)
library(survival)
library(lubridate)
library(dplyr)
library(MASS)
library(tidyr)
library(DataExplorer)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggpubr)

Final_Predict_data <- read.csv(file = '/Users/lewiswardale/Desktop/DESNT_for_Lewis/FINAL_IMPROVED_FROM_DAN/Final_Predict_Results_Dan.csv')


Predict_DESNT <- Final_Predict_data %>% 
  mutate(LPD_class = as.character(LPD_class)) %>%
  mutate(LPD_class = case_when(
    LPD_class == "LPD1" ~ "NOT_DESNT",
    LPD_class == "LPD2" ~ "NOT_DESNT",
    LPD_class == "LPD3" ~ "NOT_DESNT",
    LPD_class == "LPD4" ~ "NOT_DESNT",
    LPD_class == "LPD5" ~ "NOT_DESNT",
    LPD_class == "LPD6" ~ "NOT_DESNT",
    LPD_class == "DESNT" ~ "DESNT",
    LPD_class == "LPD8" ~ "NOT_DESNT")) %>%
  mutate(Recurrence_Event = as.character(Recurrence_Event)) %>%
  mutate(Recurrence_Event_No = case_when(
    Recurrence_Event == "FALSE" ~ "0",
    Recurrence_Event == "TRUE" ~ "1")) %>%
  mutate(Recurrence_Event_No = as.numeric(Recurrence_Event_No)) %>%
  mutate(LPD_class = as.factor(LPD_class)) %>%
  mutate(Gradegroup = as.factor(Gradegroup)) %>%
  mutate(T_Stage = as.factor(T_Stage)) %>%
  mutate(Ten_Diff = Ten_Survival_Ex_PC - Ten_Survival_ICM) %>%
  mutate(Fifteen_Diff = Fifteen_Survival_Ex_PC - Fifteen_Survival_ICM) %>%
  mutate(Ten_Diff_class = case_when(
    Ten_Diff >= 5.842 & LPD_class == "DESNT" ~ "Upper_Predict_DESNT", 
    Ten_Diff >= 5.842 & LPD_class == "NOT_DESNT" ~ "Upper_Predict_NOT_DESNT",
    Ten_Diff  < 5.842 & LPD_class == "DESNT" ~ "Lower_Predict_DESNT",
    Ten_Diff  < 5.842 & LPD_class == "NOT_DESNT" ~ "Lower_Predict_NOT_DESNT")) %>%
  mutate(Ten_Diff_class = as.factor(Ten_Diff_class)) %>%
  as.data.frame() %>%
  mutate(PSA = as.numeric(PSA)) %>%
  mutate(Age = as.numeric(Age)) %>%
  mutate(Ten_Survival_ICM = as.numeric(Ten_Survival_ICM)) %>%
  mutate(Fifteen_Survival_ICM = as.numeric(Fifteen_Survival_ICM)) %>%
  mutate(Recurrence_Time = as.numeric(Recurrence_Time)) %>%
  mutate(LPD_class = as.factor(LPD_class))

#Box Plot of DESNT status vs Survival Outcome 

Plot1 <- ggplot(Predict_DESNT, aes(x = LPD_class, y = Ten_Survival_ICM, color = LPD_class)) + 
  geom_boxplot() +
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against DESNT status") +
  labs(y = "Percentage Survival from Initial Conservative Management", x = "DESNT status", color = "DESNT status") 

Plot1 + geom_jitter(shape=16, position=position_jitter(0.2)) + theme(legend.position="none")

#Box Plot of DESNT status vs Survival Outcome 

Melted_data <- melt(Predict_DESNT, id.vars = c("Patient","Age", "Gradegroup", "PSA", "T_Stage", "Recurrence_Time", "Recurrence_Event", "LPD_class", "Dataset"), measure.vars = "Ten_Survival_ICM")

Melted_data$value = as.numeric(Melted_data$value)
Melted_data$Age = as.numeric(Melted_data$Age)
Melted_data$Recurrence_Time = as.numeric(Melted_data$Recurrence_Time)

Plot3 <- ggplot(Melted_data, aes(x = LPD_class, y = value, fill = variable) + 
  geom_boxplot(outlier.colour="red") +
  ggtitle("Survival Outcome of Initial Conservative Management against DESNT status") +
  labs(y = "Percentage Survival from Initial Conservative Management", x = "DESNT status", fill = "Survival Outcome"))

Plot3 

Man10 <- wilcox.test(Ten_Survival_ICM ~ LPD_class, data = Predict_DESNT, 
            exact = FALSE)

Man10

#Box Plot of ICM vs EXPCA 10 years 

Melted_data_10 <- melt(Predict_DESNT, id.vars = c("Patient","Age", "Gradegroup", "PSA", "T_Stage", "Recurrence_Time", "Recurrence_Event", "LPD_class", "Dataset"), measure.vars = c("Ten_Survival_ICM", "Ten_Survival_Ex_PC"))

Plot1 <- ggplot(Melted_data_10, aes(x = variable, y = value, color = variable)) + geom_boxplot() + ggtitle("Survival Outcome of Initial Conservative Management against Survival excluding Prostate Cancer deaths") + labs(y = "Percentage Survival", x = "Predict Outcome", color = "Predict Outcome")

Plot1 + geom_jitter(shape=16, position=position_jitter(0.2)) + theme(legend.position="none")

Man <- wilcox.test(value ~ variable, data = Melted_data_10, 
                     exact = FALSE)

Man


#SCATTER PLOT FOR DESNT GAMMA VS ICM SURVIVAL 

The_plot <- ggplot(Predict_DESNT, aes(x = DESNT_Gamma, y = Ten_Survival_ICM)) + geom_point() + 
  ggtitle("Survival Outcome of Initial Conservative Management against DESNT gamma") + 
  labs(y = "10 Year Percentage Survival ICM", x = "DESNT Gamma") + geom_smooth(method = "lm", se = FALSE)

cor.test(Predict_DESNT$DESNT_Gamma, Predict_DESNT$Ten_Survival_ICM, method = "spearman")

The_plot



