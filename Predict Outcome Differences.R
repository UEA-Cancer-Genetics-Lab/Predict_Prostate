
#Difference in Predict Outcome analysis

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
  mutate(Recurrence_Time = as.numeric(Recurrence_Time))

#Box Plot of DESNT status vs Survival Outcome 

Plot1 <- ggplot(Predict_DESNT, aes(x = LPD_class, y = Ten_Diff, color = LPD_class)) + 
  geom_boxplot() +
  ggtitle("10 Year Survival Outcome difference in Predict outcome against DESNT status") +
  labs(y = "Percentage difference in Survival", x = "DESNT status", color = "DESNT status") 

Plot1 + geom_jitter(shape=16, position=position_jitter(0.2)) + theme(legend.position="none")

Man10 <- wilcox.test(Ten_Diff ~ LPD_class, data = Predict_DESNT, 
                     exact = FALSE)
Man10


#Survival Analysis 

quantile(Predict_DESNT$Ten_Diff)
#0%       25%       50%       75%      100% 
#2.530271  4.448975  5.841556  7.687913 27.415518 


Predict_DESNT_Quartiles <- Predict_DESNT %>% 
  mutate(Ten_Diff_Q = case_when(
    Ten_Diff >= 4.448975 & Ten_Diff < 5.841556 ~ "Lower_Quartile10", 
    Ten_Diff >= 5.841556 & Ten_Diff <= 7.687913 ~ "Upper_Quartile10", 
    Ten_Diff > 7.687913 ~ "Above_Upper_Q10",
    Ten_Diff < 4.448975 ~ "Below_Lower_Q10"))

# COMBINED DESNT vs NOT DESNT 10

Melted_data_try <- melt(Predict_DESNT_Quartiles, id.vars = c("Patient","Age", "Gradegroup", "PSA", "T_Stage", "Recurrence_Time", "Recurrence_Event", "LPD_class", "Dataset"), measure.vars = "Ten_Diff_Q")

Melted_data_try_again <- Melted_data_try %>%
  mutate(LPD_class_Q = case_when(
    LPD_class == "NOT_DESNT" & value == "Above_Upper_Q10" ~ "ND_AUQ10",
    LPD_class == "DESNT" & value == "Above_Upper_Q10" ~ "D_AUQ10",
    LPD_class == "NOT_DESNT" & value == "Lower_Quartile10" ~ "ND_LQ10",
    LPD_class == "DESNT" & value == "Lower_Quartile10" ~ "D_LQ10",
    LPD_class == "NOT_DESNT" & value == "Upper_Quartile10" ~ "ND_UQ10",
    LPD_class == "DESNT" & value == "Upper_Quartile10" ~ "D_UQ10",
    LPD_class == "NOT_DESNT" & value == "Below_Lower_Q10" ~ "ND_BLQ10",
    LPD_class == "DESNT" & value == "Below_Lower_Q10" ~ "D_BLQ10")) %>%
  mutate(Recurrence_Event = as.character(Recurrence_Event)) %>%
  mutate(Recurrence_Event_No = case_when(
    Recurrence_Event == "FALSE" ~ "0",
    Recurrence_Event == "TRUE" ~ "1")) %>%
  mutate(Recurrence_Event_No = as.numeric(Recurrence_Event_No)) 

ICM_ALL <- survfit(Surv(Recurrence_Time, Recurrence_Event_No) ~ LPD_class_Q, data = Melted_data_try_again)
ggsurvplot(ICM_ALL, palette = c("red", "blue", "yellow", "purple", "red3", "blue3", "yellow3", "purple3"), data = Melted_data_try_again, pval = TRUE,
           risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           
           # Change ggplot2 theme 
)

surv_diff_ICM_ <- survdiff(Surv(Recurrence_Time, Recurrence_Event_No) ~ LPD_class_Q, data = Melted_data_try_again)

surv_diff_ICM_


# Scatter plot for DESNT gamma vs Differences in predict outcome

The_plot <- ggplot(Predict_DESNT, aes(x = DESNT_Gamma, y = Ten_Diff)) + geom_point() + 
  ggtitle("10 Year Survival Outcome difference in Predict outcome against DESNT gamma") + 
  labs(y = "10 Year Predict outcome difference (%)", x = "DESNT Gamma") + geom_smooth(method = "lm", se = FALSE)

cor.test(Predict_DESNT$DESNT_Gamma, Predict_DESNT$Ten_Diff, method = "spearman")

The_plot

#Scatter plot for differences vs Age

Plot1 <- ggplot(Predict_DESNT, aes(x = PSA, y = Ten_Diff, color = LPD_class)) + 
  geom_point () +
  ggtitle("10 Year Predict survival outcome differences against DESNT status") +
  labs(y = "10 Year Predict survival difference", x = "PSA", color = "DESNT status") + geom_smooth(method = "lm", se = FALSE)

Plot1

cor.test(Predict_DESNT$Ten_Diff, Predict_DESNT$PSA, method = "spearman")



