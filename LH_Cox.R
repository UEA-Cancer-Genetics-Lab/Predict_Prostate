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
library(lmtest)


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


#LH ratio test upon 2 Cox models
#DESNT + PREDICT + CLINICAL VARIABLES vs PREDICT + CLINICAL VARIABLES

DESNT_PREDICT_Variables <- coxph(Surv(Recurrence_Time, Recurrence_Event_No) ~ Age + PSA + T_Stage + Ten_Survival_ICM + LPD_class, data =  Predict_DESNT)

PREDICT_ONLY_Variables <- coxph(Surv(Recurrence_Time, Recurrence_Event_No) ~ Age + PSA + T_Stage + Ten_Survival_ICM, data =  Predict_DESNT)

lrtest(PREDICT_ONLY, DESNT_PREDICT)

#DESNT + PREDICT  vs PREDICT 

DESNT_PREDICT <- coxph(Surv(Recurrence_Time, Recurrence_Event_No) ~ Ten_Survival_ICM + LPD_class_No, data =  Predict_DESNT_No)

PREDICT_ONLY <- coxph(Surv(Recurrence_Time, Recurrence_Event_No) ~ Ten_Survival_ICM, data =  Predict_DESNT_No)

lrtest(PREDICT_ONLY, DESNT_PREDICT)

#DESNT + PREDICT  vs DESNT

DESNT_PREDICT <- coxph(Surv(Recurrence_Time, Recurrence_Event_No) ~ Ten_Survival_ICM + LPD_class_No, data =  Predict_DESNT_No)

DESNT_ONLY <- coxph(Surv(Recurrence_Time, Recurrence_Event_No) ~ LPD_class_No, data =  Predict_DESNT_No)

lrtest(DESNT_ONLY, DESNT_PREDICT)

