library(readr)
library(tidyverse)
library(prostateCancerCamcap)
library(survminer)
library(survival)
library(lubridate)
library(dplyr)
library(MASS)
library(tidyr)
library(DataExplorer)
library(magrittr)
library(ggplot2)

extra <- readRDS("DESNT2_nomogram_all_clin_data.rds")
info <- read_delim("clinical_data_with_survival_unique_patients.tsv", delim="\t")
cancermap <- read_delim("../CancerMap/clinical_data_with_survival_unique_patients.tsv", delim="\t")
camcap_int <- read_delim("../CamCap/clinical_data.tsv", delim="\t")


#Gleason score
setdiff(info$ID, extra$id)
extra_sub <- extra %>% select(c("id","Gleason_score")) 
info2 <- info %>% left_join(extra_sub,by=c("ID" = "id"))

missing_GS <- info2$ID[is.na(info2$Gleason_score)]
cm_GS <- cancermap[cancermap$ID %in% info2$ID[is.na(info2$Gleason_score)],] %>% select(c(ID,Gleason_Score))

info2$Gleason_score <- as.character(info2$Gleason_score)
info3 <- info2 %>% left_join(cm_GS, by="ID") %>% mutate(Gleason_score = case_when(is.na(Gleason_score) ~ Gleason_Score, 
                                                                         TRUE ~ Gleason_score)) %>% 
  select(-Gleason_Score)

## Can't get detailed Gleason for Glinsky

# Add Gleason grade group
info3 <- info3 %>% mutate(Gleason_Grade = case_when(
  Gleason <= 6 ~ 1,
  Gleason_score == "3+4" ~2,
  Gleason_score == "4+3" ~3,
  Gleason == 8 ~4,
  Gleason >= 9 ~5
))

# Check no ages for stockholm camcap samples - None
camcapp<- pData(camcap)
camcapp$ID <- row.names(camcapp)

cc_age <- camcap_int[camcap_int$ID %in% info2$ID[is.na(info2$Age)],] %>% 
  select(c(ID,Age.at.diag)) %>% filter(!is.na(Age.at.diag))

write_delim(info3, "improved_clindata_for_Lewis.tsv", delim="\t")

#Filtering

info_final <- info3 %>% filter(!is.na(Stage) & !is.na(Gleason_Grade) &
                   !is.na(Age)) %>% distinct(Patient_ID, .keep_all= TRUE)
dim(info_final) #363

#Further Filtering 

Final <- read.table(file = '/Users/lewiswardale/Desktop/improved_clindata_for_Lewis.tsv', sep = '\t', header = TRUE)

Final_Filtered <- Final %>%
  drop_columns(
    c("LPD1_gamma",
      "LPD2_gamma",
      "LPD3_gamma",
      "LPD4_gamma",
      "LPD5_gamma",
      "LPD6_gamma",
      "LPD8_gamma",
      "Gleason_1")) %>%
  mutate(Stage = as.character(Stage)) %>%
  mutate(Stage = case_when(
    Stage == "T1c" ~ "1",
    Stage == "T2a" ~ "2",
    Stage == "T2b" ~ "2",
    Stage == "T2c" ~ "2",
    Stage == "T3a" ~ "3",
    Stage == "T1" ~ "1",
    Stage == "T2" ~ "2",
    Stage == "T3" ~ "3",
    Stage == "T4" ~ "4",
    Stage == "T2A" ~ "2",
    Stage == "T2B" ~ "2",
    Stage == "T2C" ~ "2",
    Stage == "T3A" ~ "3",
    Stage == "T3B" ~ "3",
    Stage == "T3C" ~ "3",
    TRUE ~ Stage
  )) %>%
  dplyr::filter(!Stage == "") %>%
  dplyr::filter(!is.na(Gleason)) %>%
  dplyr::filter(!Gleason == "4") %>%
  dplyr::filter(!Age == "") %>%
  dplyr::filter(!is.na(PSA)) %>%
  filter(!duplicated(Patient_ID)) %>%
  dplyr::filter(!Patient_ID == "PR74A") %>%
  dplyr::filter(!Gleason_Grade == "") %>%
  mutate(lpd_class = as.character(lpd_class)) %>%
  mutate(lpd_class = case_when(
    lpd_class == "LPD1" ~ "LPD1",
    lpd_class == "LPD2" ~ "LPD2",
    lpd_class == "LPD3" ~ "LPD3",
    lpd_class == "LPD4" ~ "LPD4",
    lpd_class == "LPD5" ~ "LPD5",
    lpd_class == "LPD6" ~ "LPD6",
    lpd_class == "LPD7" ~ "DESNT",
    lpd_class == "LPD8" ~ "LPD8")) %>%
  mutate(dataset = as.character(dataset)) 

summary(Final_Filtered)

# Removed irrelevant columns
# Removed samples with variables missing in key clinical variables required by Predict. 

# Attempt to run Predict on each sample ----------------------------------------------------------


PredictResult <- tibble()

for(selectedPatient in unique(Final_Filtered$Patient_ID)){
  
  # Filters for each patient
  selectedPatient_Full <- Final_Filtered %>%
    filter(Patient_ID == selectedPatient) 
  
  print(selectedPatient)
  
  # Calculates Predict
  result <- predict_function(i = 15, 
                             selectedPatient_Full$Age, 
                             selectedPatient_Full$Gleason_Grade,
                             selectedPatient_Full$PSA,
                             selectedPatient_Full$Stage,
                             charlson_comorbidity = 0,
                             primaryRx = 0,
                             biopsy50 = 0)
  # Stores Predict Result
  tempVector <- c(Patient = selectedPatient, 
                  Ten_Survival_ICM = 100 - result[10,30], 
                  Fifteen_Survival_ICM = 100 - result[15,30], 
                  Ten_Survival_Ex_PC = 100 * result[10,19],
                  Fifteen_Survival_Ex_PC = 100 * result[15,19],
                  Age = selectedPatient_Full$Age,
                  Gradegroup = selectedPatient_Full$Gleason_Grade,
                  PSA = selectedPatient_Full$PSA,
                  T_Stage = selectedPatient_Full$Stage,
                  Recurrence_Time = selectedPatient_Full$recurrence_time,
                  Recurrence_Event = selectedPatient_Full$recurrence_event,
                  LPD_class = selectedPatient_Full$lpd_class,
                  Dataset = selectedPatient_Full$dataset,
                  Gleason_Score = selectedPatient_Full$Gleason_score,
                  DESNT_Gamma = selectedPatient_Full$LPD7_gamma)
  
  PredictResult <- bind_rows(PredictResult, tempVector)
  
}

write.csv(PredictResult,'/Users/lewiswardale/Desktop/Final_Predict_Results_Dan.csv', row.names =  TRUE)


# Summary of final filtered data and differences in Predict outcome....

Final_Predict_data <- read.csv(file = '/Users/lewiswardale/Desktop/DESNT_for_Lewis/FINAL_IMPROVED_FROM_DAN/Final_Predict_Results_Dan.csv')

summary(Predict_DESNT)

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
  mutate(Ten_Diff_class = as.factor(Ten_Diff_class))
