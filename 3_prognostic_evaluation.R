library(tidyverse)
library(data.table)
library(mice)
library(survival)
library(survminer)
library(cmprsk)
library(ggsurvfit)
library(tidycmprsk)
library(patchwork)
library(pROC)
rm(list = ls())
gc()
setwd(" ")

######======Prognostic analysis cohorts ======######
data <- fread("Prognostic analysis cohorts dataset.csv") 
###Variable list:
# "studyid",study_date
# "age"；"is_female"(0 male；1 female)；"region_is_urban"(0 rural；1 urban),
# "wc" (Waist Circumference, cm),"bmi"(kg/m2),
# "TG"(mmol/L),"TC"(mmol/L),"GGT"(U/L),"AST"(U/L),"ALT"(U/L),"HDLC"(mmol/L),"LDLC"(mmol/L),

###===Variables for subgroup analysis===###

##1. Alcohol-related variables:
## alcohol_category: 1 Never regular drinkers; 2 Ex-weekly drinkers; 3 Occasional or seasonal drinkers; 4 Monthly drinkers; 5 Reduced-intake drinkers; 6 Current weekly drinkers
## Only 6 current weekly drinkers were asked further questions about amount consumed for each type on a typical drinking day:
## total_alc_typ_day_g: Daily alcohol consumption (g) per occasion among weekly drinkers

## Definition of heavy alcohol consumption: (1) Ex-weekly drinkers, or (2) Reduced-intake drinkers, or (3) Current weekly drinkers who drink >15g/day for female and >30g/day for male
# Notes: In CKB, over half of ex-weekly or reduced-intake drinkers reported existing physical illness as their main reason for stopping weekly drinking
# Hence, we classified them into heavy alcohol drinkers. Subject to change according to specific cohort characteristics.
summary(as.factor(data$alcohol_CAT))

##2. diabetes diagnosis:
# fasting_glucose_x10
#diabetes diagnosis definition
# (1) History of diabetes diagnosis
# (2) Fasting glucose ≥7 mmol/L
# (3) Random glucose (fasting <8 hours) ≥11.1 mmol/L
# (4) Random glucose (fasting ≥8 hours) ≥7 mmol/L
# Notes: If criteria (3) and (4) are unavailable in the prognostic cohort, can be replaced with antidiabetic using

# In CKB
data$fbg_mmolL <- data$fasting_glucose_x10*0.1 # mmol/L
summary(data$fbg_mmolL)
data$rpg_mmolL <- data$random_glucose_x10*0.1 # mmol/L
summary(data$rpg_mmolL)
summary(as.factor(data$diabetes_diag))
data$has_diabetes_2 <- 0
data$has_diabetes_2[data$fbg_mmolL>=7] <- 1
data$has_diabetes_2[data$rpg_mmolL>=7&data$hours_since_last_ate>=8] <- 1
data$has_diabetes_2[data$rpg_mmolL>=11.1&data$hours_since_last_ate<8] <- 1
data$has_diabetes_2[data$diabetes_diag==1] <- 1
data$has_diabetes_2[data$is_diabetic==1] <- 1
summary(as.factor(data$has_diabetes_2))

###3. fatty liver subgroups: FLI >=60####
# Fatty liver definition: FLI > 60: Fatty liver；FLI ≤ 60: Non-fatty liver
# FLI calculation: TG (mg/dl), BMI (kg/m²), GGT (U/L), WC (cm)
data$tg_mgdl <- 88.57*data$TG
data$sum <- 0.953*log(data$tg_mgdl) + 0.139*data$bmi + 0.718*log(data$GGT) + 0.053*data$wc - 15.745
data$fli <- exp(data$sum)/(1+exp(data$sum))*100
data$fli_60 <- 0
data$fli_60[data$fli >= 60] <- 1

###===Diseases endpoints for prognostic analysis===###
# All-cause death, liver related diseases event and death(ICD-10: C22.0,C22.7,K70.1,K70.2,K70.3,K70.4,K70.9,K71.7,K71.9,K72.0,K72.1,K72.9,
# K73.0,K73.2,K73.8,K73.9","K74.0,K74.1,K74.2,K75.8,K75.9,K76.6,K76.7,K76.8,K76.9,I85.0,I85.9,I98.2,I98.3)
# Censoring Date：2022/12/31（CKB）

#  baseline survey date
summary(data$study_date)
# LRE (liver related event) ，
data$lre
data$lre_date
# LRM (liver related mortality)
data$lrm
data$lrm_date

##Follow-up duration was calculated from baseline survey date to the earliest occurrence of outcome event  or censoring.
var <- c("lre","lrm")
date_var <- c("lre_date","lrm_date")
for (i in 1:length(var)) {
        time_var=  as.numeric(as.Date(substr(data[[date_var[i]]], 1, 10), "%Y-%m-%d") - as.Date(substr(data$study_date, 1, 10), "%Y-%m-%d")) / 365.25
        data <- data %>%
                mutate(!!sym(paste0(var[i], "_time")) := time_var )
}


###***Table 1. Baseline characteristics of the study population***###
# Variables
all_vars <- c("age","is_female","bmi","has_diabetes","AST","ALT","GGT","TC","TG","LDLC","HDLC")
cat_vars <- c("is_female","has_diabetes")
# TableOne
library(tableone)
tab1 <- CreateTableOne(vars = all_vars,data = data,factorVars = cat_vars,addOverall = TRUE)
tab1_output<- print(tab1, showAllLevels = TRUE,nonnormal = NULL,digits = 1,catDigits = 1,  formatOptions = list(big.mark = ","))
tab1_output <- as.data.frame.matrix(tab1_output)
write.csv(tab1_output,file="Table 1.csv")



#####===== CKB cirrhosis risk score  CCS=====#####
## the CKB Cirrhosis Score (CCS) using a penalized logistic regression
## Predictors: age (years), diabetes（0：No；1：Yes）, AST (U/L), GGT (spline variable with 3 knots, U/L), TC（mmol/L）, HDL-C（mmol/L）
#  and the interaction term between sex （0：Men；1：Women）and BMI (spline variable with 3 knots, Figure 2).
final_model_CCS <- function(age = NA,is_female = NA,
                            has_diabetes = NA,bmi = NA,tc = NA,hdl = NA,ast = NA,ggt = NA){
        final_model_score = -5.239626+0.046680734*age-5.6458036*is_female+0.46411484*has_diabetes-0.16098654* bmi+0.0026687438*pmax(bmi-20.200001,0)^3-0.0049893916*pmax(bmi-24.200001,0)^3+0.0023206478*pmax(bmi-28.799999,0)^3+
                is_female*(0.26032348* bmi-0.0031969233*pmax(bmi-20.200001,0)^3+0.0059768578*pmax(bmi-24.200001,0)^3-0.0027799345*pmax(bmi-28.799999,0)^3)-0.1056094*tc -0.40231523*hdl  +0.013602785*ast+0.10879714* ggt-0.0001017148*pmax(ggt-13.4,0)^3+0.00013169036*pmax(ggt-22.3,0)^3-2.9975555e-05*pmax(ggt-52.5,0)^3
        final_model_prop = 1/(1+exp(-final_model_score))
        return(final_model_prop)
}
data$CCS_exp_risk <- final_model_CCS(
        age= data$age,
        is_female = as.numeric(data$is_female),
        bmi=data$bmi,
        ggt = data$GGT,
        ast=data$AST,
        tc=data$TC,
        hdl = data$HDLC,
        has_diabetes = as.numeric(data$has_diabetes)
)
##Risk scores cut-off:  2%, 10%
data <- data %>% mutate(CCS_exp_risk_cat = cut(CCS_exp_risk, breaks = c(0, 0.02, 0.1, 1), labels = c("Low", "Medium","High"), include.lowest = TRUE)%>%factor(levels = c("Low", "Medium","High")))
#< 2%, Low risk group; 2-10% Median risk group；>=10% High risk group
table(data$CCS_exp_risk_cat)
data$CCS_exp_risk_z <- scale(data$CCS_exp_risk, scale = TRUE, center = TRUE)


###***Table 2. Prognostic performance of CCS ***###
##Definition of competing risk outcomes
##All-cause mortality
data$death
##All-cause mortality follow-up time was calculated from baseline survey to censoring
data$time_all <- as.numeric(as.Date(substr(data$censoring_date2, 1, 10), "%Y-%m-%d") - as.Date(substr(data$study_date, 1, 10), "%Y-%m-%d")) / 365.25
##Liver related mortality: lrm
summary(as.factor(data$lrm)) 

### Competing risk outcomes：event=1: censored, event=2: outcome of interest, event=3: competing event
### Competing risk follow-up time ：If the event of interest occurs, the time of the event of interest is taken; otherwise, the time of all-cause death is taken
var <- c("lre","lrm")
for ( i in 1:length(var)){
        ##Define the endpoints with competing outcomes
        data$lr_status <- 1
        data$lr_status[data$death == 1 & data$lrm == 0] <- 3
        data$lr_status[data[[var[i]]]==1] <- 2
        data$lr_status_time <- ifelse(data$lr_status==2, data[[paste0(var[i],"_time")]],data$time_all)
        data[[paste0(var[i],"_status_time")]] <- data$lr_status_time
        data[[paste0(var[i],"_status")]] <- as.factor(data$lr_status)
}

###=====SHR: Fine-Gray model, considering the competing risk outcomes=======###
library(cmprsk)
library(tidyr)
competing_ends <- c("lre_status","lrm_status")
HR_FG <- NA
for ( end in competing_ends){
        data$lr_status <-  as.factor(data[[end]])
        data$lr_status_time <- data[[paste0(end,"_time")]]
        
        crr_lr <- crr(Surv(lr_status_time, lr_status) ~ CCS_exp_risk_cat, data = data, cencode = 1,failcode=2)
        crr_lr_cat <- tidy(crr_lr) %>% as.data.frame %>%
                mutate(`HR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", exp(estimate), exp(estimate-1.96*std.error), exp(estimate+1.96*std.error)))
        HR <- data.frame("HR_FG"=crr_lr_cat$'HR (95% CI)',"Var"=end)
        HR_FG <- rbind(HR_FG,HR)
}
HR_FG 

###=====Wolber’s C: Fine-Gray  model=======###
### Wolber's C-index ###
library(pec)
library(riskRegression)
competing_ends <- c("lre_status","lrm_status")
Cindexes <- NA
for ( end in competing_ends){
        data$lr_status <-  data[[end]]
        data$lr_status_time <- data[[paste0(end,"_time")]]
        ##Define the endpoints with competing outcomes
        fgr_lr <- FGR(Hist(lr_status_time, lr_status) ~ CCS_exp_risk_z, data=data, cause=2)
        cindex_result<- cindex(list("Fine-gray model" = fgr_lr),
                               formula = Hist(lr_status_time,lr_status) ~ 1,
                               cens.model="marginal",
                               data=data, cause=2,confInt = TRUE)       
        Cindex <- do.call(cbind,cindex_result$AppCindex) %>% data.frame()
        Cindex$var <- end
        Cindexes <- rbind(Cindexes,Cindex)
}

Cindexes

##Time dependent -AUROC considered competing risk 
library("timeROC")
competing_ends <- c("lre_status","lrm_status")
AUC_table_competing <-NA
for (  end in competing_ends){
        data$lr_status <-  data[[end]]
        data$lr_status_time <- data[[paste0(end,"_time")]]
        select_var <-c("studyid","lr_status","lr_status_time","CCS_exp_risk_z")
        data1 <-data %>% select(all_of(select_var))
        
        roc1<-timeROC(T=data1$lr_status_time,
                      delta=data1$lr_status, marker=data1$CCS_exp_risk_z,
                      cause=2, weighting="marginal",
                      times=c(2,5,10),
                      ROC = TRUE, iid =F)
        
        ##AUC table###
        tabler1 <- data.frame(c(sprintf("%.3f",roc1$AUC_2)))
        rownames(tabler1) <-paste0(c(2,5,10)[1:nrow(tabler1)]," years",end)
        AUC_table_competing <-rbind(AUC_table_competing,tabler1)
}
write.csv(AUC_table_competing,file="Table 2.AUC_table_competing.csv")



###***Figure 4. Cumulative incidence of LRE and LRM in the CKB prognostic cohort and the cohort prognostic cohort***###
summary(data$lre_status);summary(data$lrm_status) #event=1: censored, event=2: outcome of interest, event=3: competing event
summary(data$lre_status_time) ;summary(data$lrm_status_time)
endpoints <- c("lre_status","lrm_status")
##CCS risk category 
summary(data$CCS_exp_risk_cat)

my_colors <- c("#51AD43", "#FB7B00",  "#DF0906")
### CIF: Liver related event LRE ###
cif_lre <- cuminc(Surv(lrm_status_time, lre_status) ~ CCS_exp_risk_cat, data = data,cencode=1) %>%
        ggcuminc(outcome = "2",
                 theme = theme_classic(base_size = 12),
                 size = 1) +
        scale_x_continuous(
                name = "Time (years)", 
                breaks = seq(0, 20, by = 5),  
                limits = c(0, 20)
        ) +
        scale_y_continuous(
                name = "Cumulative Incidence (%)",
                breaks = seq(0, 0.08, by = 0.02), 
                labels = sprintf("%.1f",seq(0, 0.08, by = 0.02)*100),
                limits = c(0, 0.08)
        )+         
        theme(legend.position = "inside", 
              legend.position.inside = c(0.25, 0.9),
              legend.key.size = unit(1.5, "lines"),
              text = element_text(color = "black"),  
              axis.text = element_text(color = "black"), 
              axis.title = element_text(color = "black"), 
              legend.text = element_text(color = "black") 
        )+ 
        labs(x = "Times (years)",
             title = "")+
        scale_color_manual(values = my_colors, limits = levels(data$CCS_exp_risk_cat)) +
        scale_fill_manual(values = my_colors, limits = levels(data$CCS_exp_risk_cat),
                          labels = c("Low risk", "Medium risk", "High risk"))

### CIF : LRM considered competing risks
cif_lrm <- cuminc(Surv(lrm_time, lrm_status) ~ CCS_exp_risk_cat, data = data,cencode=1) %>%
        ggcuminc(outcome = "2",
                 theme = theme_classic(base_size = 12),
                 size = 1) +
        scale_x_continuous(
                name = "Time (years)",  
                breaks = seq(0, 20, by = 5), 
                limits = c(0, 20)
        ) +
        scale_y_continuous(
                name = "Cumulative Incidence (%)",
                breaks = seq(0, 0.02, by = 0.01), 
                labels = sprintf("%.1f",seq(0, 0.02, by = 0.01)*100),
                limits = c(0, 0.02)#,
        )+         
        theme(legend.position = "inside", 
              legend.position.inside = c(0.25, 0.9),
              legend.key.size = unit(1.5, "lines"),
              text = element_text(color = "black"),  
              axis.text = element_text(color = "black"), 
              axis.title = element_text(color = "black"), 
              legend.text = element_text(color = "black") 
        )+ 
        labs(x = "Times (years)",
             title = "")+
        scale_color_manual(values = my_colors, limits = levels(data$CCS_exp_risk_cat)) +
        scale_fill_manual(values = my_colors, limits = levels(data$CCS_exp_risk_cat),
                          labels = c("Low risk", "Medium risk", "High risk"))

plots_end <- ggarrange(
        labels = c("(A) LRE","(B) LRM "),
        plotlist = list(cif_lre,  cif_lrm),
        nrow = 2, ncol =1, align = "hv", common.legend = TRUE, legend = "bottom"
)
ggsave(filename = paste0("Output/Figure 4.CIF_LRE_LRM.pdf"), 
       plot = plots_end, width = 5.5, height = 11, units = "in", dpi = 300)

