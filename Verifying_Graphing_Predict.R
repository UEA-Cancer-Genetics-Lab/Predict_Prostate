#Graphing Predict results

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

# Melted data set 

Melted_data <- Predict_DESNT %>%
  melt(Predict_DESNT, id.vars = c("Patient","Age", "Gradegroup", "PSA", "T_Stage", "Recurrence_Time", "Recurrence_Event", "LPD_class", "Dataset"), measure.vars = "Ten_Survival_ICM") %>%
  mutate(value = as.factor(value)) %>%
  mutate(Age = as.factor(Age)) %>%
  mutate(Recurrence_Time = as.factor(Recurrence_Time)) 


# ICM 10 Year vs AGE 

ICM_Ten_Age_ALL <- ggplot(data = Predict_DESNT, aes(x = Age , y = Ten_Survival_ICM)) + geom_point() + theme_minimal() + geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against Age at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management") + scale_x_continuous(limits=c(35,75) ,breaks=c(35, 45, 55, 65, 75)) + scale_y_continuous(limits=c(50,100))

ICM_Ten_Age_ALL

cor.test(Predict_DESNT$Age, Predict_DESNT$Ten_Survival_ICM, method = "spearman")

# ICM 10 Year vs PSA

ICM_Ten_PSA_ALL <- ggplot(data = Predict_DESNT, aes(x = PSA , y = Ten_Survival_ICM)) + geom_point() + theme_minimal() + geom_smooth(method = lm, se = FALSE, fullrange = TRUE) + 
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against PSA at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management")

ICM_Ten_PSA_ALL

cor.test(Predict_DESNT$PSA, Predict_DESNT$Ten_Survival_ICM, method = "spearman")

# ICM 10 Year vs PSA + LPD)

ICM_Ten_PSA_LPD <- ggplot(data = Predict_DESNT, aes(x = PSA , y = Ten_Survival_ICM, color = LPD_class)) + geom_point() + theme_minimal() +
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against PSA at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management", color = "Assigned_LPD")

ICM_Ten_PSA_LPD

# ICM 10 Year vs GradeGroup  (ALL)

ICM_Ten_GRADE_ALL <- ggplot(data = Predict_DESNT, aes(x = Age , y = Ten_Survival_ICM, color = Gradegroup)) + geom_point() + theme_minimal() + 
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against Age and Grade at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management", color = "Grade Group") + scale_x_continuous(limits=c(35,75) ,breaks=c(35, 45, 55, 65, 75))

ICM_Ten_GRADE_ALL

# ICM 10 Year vs T_Stage  (ALL)

ICM_Ten_TStage_ALL <- ggplot(data = Predict_DESNT, aes(x = Age , y = Ten_Survival_ICM, color = T_Stage)) + geom_point() + theme_minimal() + 
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against Age and T-Stage at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management", color = "T-Stage") + scale_x_continuous(limits=c(35,75) ,breaks=c(35, 45, 55, 65, 75))

ICM_Ten_TStage_ALL

# Age vs Recurrence event 

ICM_Age_ALL_RE_MELT <- ggplot(data = Melted_data, aes(x = Age , y = value, color = Recurrence_Event)) + geom_point() + theme_minimal() + 
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against Age for Recurrence Event at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management", color = "Recurrence Event Outcome") + scale_x_continuous(limits=c(35,75) ,breaks=c(35, 45, 55, 65, 75))

ICM_Age_ALL_RE_MELT

# Recurrence Time vs Survival 

ICM_Ten_RT_ALL <- ggplot(data = Predict_DESNT, aes(x = Recurrence_Time , y = Ten_Survival_ICM, color = Recurrence_Event)) + geom_point() + theme_minimal() + 
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against Age and T-Stage at Diagnosis") +
  labs(y = "Percentage Survival from Initial Conservative Management", color = "Recurrence Status") + scale_x_continuous(limits=c(0,100) ,breaks=c(0, 20,40 ,60 , 80, 100))

ICM_Ten_RT_ALL

# DESNT GAMMA vs Survival 

ICM_Ten_Gamma <- ggplot(data = Predict_DESNT, aes(x = DESNT_Gamma , y = Ten_Survival_ICM)) + geom_point() + theme_minimal() + 
  ggtitle("10 Year Survival Outcome of Initial Conservative Management against DESNT Gamma") +
  labs(y = "Percentage Survival from Initial Conservative Management") 

ICM_Ten_Gamma


