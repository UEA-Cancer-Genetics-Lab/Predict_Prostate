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
  
# Kaplan-Meier DESNT vs NOT DESNT

f1 <- survfit(Surv(Recurrence_Time, Recurrence_Event_No) ~ LPD_class, data = Predict_DESNT)

ggsurvplot(f1, data = Predict_DESNT, pval = TRUE, conf.int = TRUE,
risk.table = TRUE, # Add risk table
risk.table.col = "strata", # Change risk table color by groups
linetype = "strata", # Change line type by groups
surv.median.line = "hv", # Specify median survival
ggtheme = theme_bw(), # Change ggplot2 theme
palette = c("#E7B800", "#2E9FDF"))

surv_diff_DESNT <- survdiff(Surv(Recurrence_Time, Recurrence_Event_No) ~ LPD_class, data = Predict_DESNT)

surv_diff_DESNT

#Log_test
#Chisq= 48.6  on 1 degrees of freedom, p= 3e-12 

# Kaplan-Meier on Predict Conserved management 

quantile(Predict_DESNT$Ten_Survival_ICM)

# 0%      25%      50%      75%     100% 
# 54.41731 81.05083 87.33297 91.02882 96.53431 

quantile(Predict_DESNT$Ten_Survival_Ex_PC)

# 0%      25%      50%      75%     100% 
# 68.36567 89.18849 93.88986 96.44686 99.58791 
# 68.36567 89.28694 93.52898 96.44198 99.58791 

Predict_DESNT_Quartiles <- Predict_DESNT %>% 
  mutate(Ten_Survival_ICM_Q = case_when(
    Ten_Survival_ICM >= 81.05083 & Ten_Survival_ICM < 87.33297 ~ "Lower_Quartile10", 
    Ten_Survival_ICM >= 87.33297 & Ten_Survival_ICM <= 91.02882 ~ "Upper_Quartile10", 
    Ten_Survival_ICM > 91.02882 ~ "Above_Upper_Q10",
    Ten_Survival_ICM < 81.05083 ~ "Below_Lower_Q10")) %>%
  mutate(Ten_Survival_Ex_PC_Q = case_when(
    Ten_Survival_Ex_PC >= 89.28694 & Ten_Survival_Ex_PC < 93.52898 ~ "Lower_Quartile10_EXPC", 
    Ten_Survival_Ex_PC >= 93.52898 & Ten_Survival_Ex_PC <= 96.44198 ~ "Upper_Quartile10_EXPC", 
    Ten_Survival_Ex_PC > 96.44198 ~ "Above_Upper_Q_EXPC10",
    Ten_Survival_Ex_PC < 89.28694 ~ "Below_Lower_Q_EXPC10"))


f3 <- survfit(Surv(Recurrence_Time, Recurrence_Event_No) ~ Ten_Survival_ICM_Q, data = Predict_DESNT_Quartiles)

ggsurvplot(f3, data = Predict_DESNT_Quartiles, pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme 
)

surv_diff_QICM <- survdiff(Surv(Recurrence_Time, Recurrence_Event_No) ~ Ten_Survival_ICM_Q, data = Predict_DESNT_Quartiles)

surv_diff_QICM

#Log_test_10yr_ICM
#Chisq= 31  on 3 degrees of freedom, p= 9e-07 

# COMBINED DESNT vs NOT DESNT 

Melted_data_try <- melt(Predict_DESNT_Quartiles, id.vars = c("Patient","Age", "Gradegroup", "PSA", "T_Stage", "Recurrence_Time", "Recurrence_Event", "LPD_class", "Dataset"), measure.vars = c("Ten_Survival_ICM_Q"))

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
    


ICM_vs_ExPc <- survfit(Surv(Recurrence_Time, Recurrence_Event_No) ~ Survival_Type_Quartiles , data = Melted_data_ICM_vs_EXPC)

ggsurvplot(ICM_vs_ExPc, palette = c("red", "blue", "yellow", "purple", "red1", "blue1", "yellow1", "purple1" ), data = Melted_data_ICM_vs_EXPC, pval = TRUE,
           risk.table = TRUE
)

# COMBINED 10 year ICM vs 10 year Excluding PCA deaths  


Melted_data_ICM_vs_EXPC_10 <- melt(Predict_DESNT_Quartiles, id.vars = c("Patient","Age", "Gradegroup", "PSA", "T_Stage", "Recurrence_Time", "Recurrence_Event", "LPD_class", "Dataset", "DESNT_Gamma"), measure.vars = c("Ten_Survival_ICM_Q", "Ten_Survival_Ex_PC_Q")) %>%
  mutate(Recurrence_Event = as.character(Recurrence_Event)) %>%
  mutate(Recurrence_Event_No = case_when(
    Recurrence_Event == "FALSE" ~ "0",
    Recurrence_Event == "TRUE" ~ "1")) %>%
  mutate(Recurrence_Event_No = as.numeric(Recurrence_Event_No)) %>%
  mutate(Survival_Type_Quartiles = case_when(
    variable == "Ten_Survival_Ex_PC_Q" & value == "Above_Upper_Q_EXPC10" ~ "EXPCA_AUQ10",
    variable == "Ten_Survival_Ex_PC_Q" & value == "Below_Lower_Q_EXPC10" ~ "EXPCA_BLQ10",
    variable == "Ten_Survival_Ex_PC_Q" & value == "Lower_Quartile10_EXPC" ~ "EXPCA_LQ10",
    variable == "Ten_Survival_Ex_PC_Q" & value == "Upper_Quartile10_EXPC" ~ "EXPCA_UQ10",
  ))

ICM_vs_ExPc_10 <- survfit(Surv(Recurrence_Time, Recurrence_Event_No) ~ Survival_Type_Quartiles , data = Melted_data_ICM_vs_EXPC_10)

ggsurvplot(ICM_vs_ExPc_10, palette = c("red", "blue", "yellow", "purple", "red1", "blue1", "yellow1", "purple1" ), data = Melted_data_ICM_vs_EXPC_10, pval = TRUE,
           risk.table = TRUE
)

surv_diff_ICM_ <- survdiff(Surv(Recurrence_Time, Recurrence_Event_No) ~ Survival_Type_Quartiles, data =  Melted_data_ICM_vs_EXPC_10)

surv_diff_ICM_
