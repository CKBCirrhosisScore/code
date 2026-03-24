rm(list = ls())
gc()

library(tidyverse)
library(survey)
library(foreign)
library(boot)
library(WeightedROC)
library(pROC)
library(rms)
library(gam)
library(xgboost)
library(tableone)
library(PredictABEL)
library(nricens)
library(scales)
setwd("")

#### ============================ Load data set and clean variables ========================= ####
load("NHANES data 2017-2020.RData")

# nhanes0: 4476
nhanes0 <- nhanes0[!is.na(nhanes0$WTSAFPRP) & nhanes0$WTSAFPRP > 0,] %>%
        mutate(inAnalysis = (age >= 18 & is_pregnant != 1 & vcte_valid == 1 & !is.na(vcte_valid) &
                                     !is.na(median_stiff) & !is.na(median_fat) &
                                     !is.na(age) & !is.na(is_female) & !is.na(bmi) & !is.na(wc) & !is.na(hip) & 
                                     !is.na(ast) & !is.na(alt) & !is.na(ggt) & !is.na(alb) & !is.na(alp) & !is.na(tbil) & !is.na(sodium) & !is.na(plt) & 
                                     !is.na(tc) & !is.na(tg) & !is.na(hdl) & !is.na(ldl) &
                                     !is.na(diabetes_diag) & !is.na(fbg) & !is.na(hypertension_diag) & !is.na(sbp) & !is.na(dbp)))

summary(as.factor(nhanes0$diabetes_diag)) 
nhanes0$has_diabetes <- 0
nhanes0$has_diabetes[nhanes0$diabetes_diag == 1] <- 1
nhanes0$has_diabetes[nhanes0$fbg >= 7] <- 1
summary(as.factor(nhanes0$has_diabetes))
nhanes0$fbg_mgdl <- nhanes0$fbg * 18

nhanes0$has_hyperglycemia <- 0
nhanes0$has_hyperglycemia[nhanes0$diabetes_diag == 1] <- 1
nhanes0$has_hyperglycemia[nhanes0$fbg >= 5.56] <- 1
summary(as.factor(nhanes0$has_hyperglycemia))

nhanes0$has_htn <- 0
nhanes0$has_htn[nhanes0$hypertension_diag == 1] <- 1
nhanes0$has_htn[nhanes0$sbp >= 130] <- 1
nhanes0$has_htn[nhanes0$dbp >= 80] <- 1
summary(as.factor(nhanes0$has_htn))

summary(nhanes0$avg_alco)
nhanes0$heavy_alco <- 0
nhanes0$heavy_alco[(nhanes0$avg_alco >= 2 & nhanes0$is_female == 0) | (nhanes0$avg_alco >= 1 & nhanes0$is_female == 1)] <- 1
summary(as.factor(nhanes0$heavy_alco))

nhanes0$fibrosis_f4_ft <- NA
nhanes0$fibrosis_f4_ft[!is.na(nhanes0$median_stiff)] <- 0
nhanes0$fibrosis_f4_ft[nhanes0$median_stiff>=13.5] <- 1
summary(as.factor(nhanes0$fibrosis_f4_ft)) 

nhanes0$fibrosis_f4_fs <- NA
nhanes0$fibrosis_f4_fs[!is.na(nhanes0$median_stiff)] <- 0
nhanes0$fibrosis_f4_fs[nhanes0$median_stiff>=15] <- 1
summary(as.factor(nhanes0$fibrosis_f4_fs))

nhanes0$fibrosis_f4_plt <- NA
nhanes0$fibrosis_f4_plt[!is.na(nhanes0$median_stiff) & !is.na(nhanes0$plt)] <- 0
nhanes0$fibrosis_f4_plt[nhanes0$median_stiff > 20 & nhanes0$plt < 150] <- 1
summary(as.factor(nhanes0$fibrosis_f4_plt))
#    0    1 NA's 
# 4217   11  248

nhanes0$region_is_urban <- 0
nhanes0$tc_z <- scale(nhanes0$tc, center = TRUE, scale = TRUE)
nhanes0$hdl_z <- scale(nhanes0$hdl, center = TRUE, scale = TRUE)
nhanes0$ast_z <- scale(nhanes0$ast, center = TRUE, scale = TRUE)
nhanes0$ggt_z <- scale(nhanes0$ggt, center = TRUE, scale = TRUE)

nhanes0 <- nhanes0 %>%
        mutate(lp_ccs = -5.239626+0.10879714* ggt-0.0001017148*pmax(ggt-13.4,0)^3+0.00013169036*pmax(ggt-22.3,0)^3-2.9975555e-05*pmax(ggt-52.5,0)^3+0.046680734*age-0.1056094*tc-5.6458036*is_female+0.013602785*ast+0.46411484*has_diabetes-0.40231523*hdl-0.16098654* bmi+0.0026687438*pmax(bmi-20.200001,0)^3-0.0049893916*pmax(bmi-24.200001,0)^3+0.0023206478*pmax(bmi-28.799999,0)^3+is_female*(0.26032348* bmi-0.0031969233*pmax(bmi-20.200001,0)^3+0.0059768578*pmax(bmi-24.200001,0)^3-0.0027799345*pmax(bmi-28.799999,0)^3),
               p_ccs = exp(lp_ccs)/(exp(lp_ccs)+ 1),
               fib4 = age*ast/(plt*sqrt(alt)),
               apri = ast*100/(40*plt),
               lp_apri = 2.411 + 0.100*log(ast)-0.436*log(plt),
               p_apri = exp(lp_apri)/(exp(lp_apri) + 1),
               forn = 7.811 - 3.131*log(plt) + 0.781*log(ggt) + 3.467*log(age)-0.014*tc,
               p_forn = exp(forn)/(exp(forn) + 1),
               nfs = -1.675 + (0.037 * age) + (0.094 * bmi) + (1.13 * has_hyperglycemia) + (0.99 * ast / alt) - (0.013 * plt) - (0.66 * alb / 10),
               p_nfs = exp(nfs)/(exp(nfs) + 1),
               ccs_cat = cut(p_ccs, breaks = c(-Inf, 0.02, 0.1, Inf), include.lowest = TRUE, right = FALSE),
               fib20 = -12.350448 + 0.014748012 * ggt - 1.0548047e-05 * pmax(ggt - 13.4, 0)^3 + 1.3656577e-05 * pmax(ggt - 22.3, 0)^3 - 3.1085304e-06 * pmax(ggt - 52.5, 0)^3 + 0.056497217 * age - 0.15712383 * tc + 0.64774771 * is_female + 0.013840679 * ast + 0.30925669 * has_diabetes - 0.26550003 * hdl + 0.086576651 * bmi + 0.00044857635 * pmax(bmi - 20.200001, 0)^3 - 0.00083864291 * pmax(bmi - 24.200001, 0)^3 + 0.00039006656 * pmax(bmi - 28.799999, 0)^3 + is_female * (-0.022050423 * bmi - 0.00076288205 * pmax(bmi - 20.200001, 0)^3 + 0.001426258 * pmax(bmi - 24.200001, 0)^3 - 0.00066337598 * pmax(bmi - 28.799999, 0)^3),
               p_fib20 = 1/(1+exp(-fib20)))

## using complete-case analysis
nhanes <- nhanes0[nhanes0$inAnalysis, ] # 3132

### exact LiverRisk score: downloaded from the website (https://www.liverriskscore.com/)
# nhanes_to_liverrisk <- nhanes %>% 
#         mutate(sex = ifelse(is_female == 1, "Female", "Male")) %>%
#         select(sex, age, fbg, tc, ast, alt, ggt, plt)
# write.csv(nhanes_to_liverrisk, "nhanes_exact_liverrisk.csv", row.names = FALSE)

nhanes_exact_liverrisk <- read.csv("nhanes_exact_liverrisk_output.csv")
summary(nhanes_exact_liverrisk)
nhanes$liverrisk <- nhanes_exact_liverrisk$score

### eSL
nhanesX <- nhanes[, c("age",  "whr", "ast", "ggt", "alt", "hdl", "tc", "tg", 
                      "fbg", "sbp", "dbp", "region_is_urban", "has_diabetes", 
                      "has_htn", "is_female", "bmi")]

names(nhanesX) <- paste0(names(nhanesX), "_3r")
names(nhanesX)[c(15, 12)] <- c("is_female", "region_is_urban")
names(nhanesX)[9] <- "rpg_3r"

load("2025-01-16_SL_8_models_F4_16pred.rData")

nhanes$eSL <- predict(SL_8_F4, nhanesX)$pred

### Add liverrisk and eSL back in nhanes0
nhanes0 <- nhanes0 %>% left_join(nhanes %>% select(SEQN, liverrisk, eSL), by = "SEQN")

## multi-stage sampling design
design = svydesign(id = ~ SDMVPSU, 
                   strata = ~SDMVSTRA, 
                   weights = ~WTSAFPRP, nest=TRUE, 
                   data = nhanes0)
nhanes_design <- subset(design, inAnalysis)

#### ============================ Table 1 ========================= ####
### Overall ###
svyCreateTableOne(vars = c("age", "is_female", "bmi", "has_diabetes", "ast", "alt", "ggt", "tc", "tg", "ldl", "hdl","median_stiff", "fibrosis_f4_ft"), 
                  strata = "fibrosis_f4_ft",
                  data = nhanes_design, 
                  factorVars = c("is_female", "has_diabetes", "fibrosis_f4_ft"))
table(nhanes$is_female)
table(nhanes$has_diabetes)

### Stratified by cirrhosis status ###
svyCreateTableOne(vars = c("age", "is_female", "bmi", "has_diabetes", "ast", "alt", "ggt", "tc", "tg", "ldl", "hdl","median_stiff", "fibrosis_f4_ft"), 
                  data = nhanes_design, 
                  factorVars = c("is_female", "has_diabetes", "fibrosis_f4_ft"))
table(nhanes$fibrosis_f4_ft)
tapply(nhanes$is_female, nhanes$fibrosis_f4_ft, table)
tapply(nhanes$has_diabetes, nhanes$fibrosis_f4_ft, table)

#### =============================== AUROC =================================== ####

### Refitted value of different models ###
# Fit predefined models
fit_liverrisk <- svyglm(fibrosis_f4_ft ~ age + is_female + fbg + tc + ast + alt + ggt + plt,
                        family = quasibinomial, design = nhanes_design)
nhanes$pred_liverrisk <- predict.glm(fit_liverrisk, data = nhanes, type = "response")
summary(nhanes$pred_liverrisk)

fit_liverrisk_lsm <- svyglm(median_stiff ~ age + is_female + fbg + tc + ast + alt + ggt + plt,
                            family = "gaussian", design = nhanes_design)
nhanes$pred_liverrisk_lsm <- predict.glm(fit_liverrisk_lsm, data = nhanes)
summary(nhanes$pred_liverrisk_lsm)

# According to preliminary analysis, this combination has the highest AUROC
fit_liverpro <- svyglm(fibrosis_f4_ft ~ age + is_female +  tc + ast + ggt + alp + alb + tbil,
                       family = quasibinomial, design = nhanes_design)
nhanes$pred_liverpro <- predict.glm(fit_liverpro, data = nhanes, type = "response")
summary(nhanes$pred_liverpro)

# Function to calculate weighted AUC for multiple predictors
calc_AUC_weighted = function(preds, response, weights) {
        sapply(preds, function(x) {
                WeightedAUC(WeightedROC(x, response, weights))
        })
}

# Function to calculate bootstrap weighted AUCs for multiple scores
boot_AUC_weighted = function(preds, response, weights, B = 1000) {
        run_boot = function(){
                boot_idx = sample(length(response), replace = TRUE)
                sapply(preds, function(x) {
                        WeightedAUC(WeightedROC(x[boot_idx], response[boot_idx], weights[boot_idx]))
                })
        }
        return(t(replicate(B, run_boot())))
}

# Set seed for reproducibility
set.seed(1)

# Compute AUCs and bootstrap distributions
AUC_nhanes = calc_AUC_weighted(nhanes[, c("p_ccs", "eSL", "apri", "fib4", "forn", "nfs", "liverrisk", "pred_liverrisk", "pred_liverrisk_lsm", "pred_liverpro")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)
boot_AUC_nhanes = boot_AUC_weighted(nhanes[, c("p_ccs", "eSL", "apri", "fib4", "forn", "nfs", "liverrisk", "pred_liverrisk", "pred_liverrisk_lsm", "pred_liverpro")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

# Compute confidence intervals
AUC_CI_df = cbind(
        "AUROC" = AUC_nhanes, 
        "AUROC 2.5%" = apply(boot_AUC_nhanes, 2, quantile, 0.025), 
        "AUROC 97.5%" = apply(boot_AUC_nhanes, 2, quantile, 0.975)
) %>% as.data.frame()
AUC_CI_df$auroc <- sprintf("%.3f (%.3f, %.3f)", AUC_CI_df$AUROC, AUC_CI_df$`AUROC 2.5%`, AUC_CI_df$`AUROC 97.5%`)
AUC_CI_df

#### ============================ Calibration ========================= ####
# CCS deciles for each validation set
p_ccs_dec = as.numeric(cut(nhanes$p_ccs, 
                           breaks = c(-Inf, as.numeric(quantile(nhanes$p_ccs, seq(0.1, 0.9, by = 0.1))), Inf),
                           labels = c(1:(length(quantile(nhanes$p_ccs, seq(0.1, 0.9, by = 0.1)))+1))))

### Brier score ###
# Function to calculate Brier score for multiple scores, 
# while accounting for sampling weights
calc_Brier_weighted = function(preds, response, weights) {
        sapply(preds, function(x) {
                Brier = weighted.mean((x - response)^2, weights)
        })
}

Brier_nhanes = calc_Brier_weighted(nhanes[, c("p_ccs", "eSL")], 
                                   nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

# Function to calculate bootstrap weighted Brier score for multiple scores
boot_Brier_weighted = function(preds, response,  weights, B = 1000) {
        run_boot = function(){
                # Sample indices with replacement
                boot_idx = sample(length(response), replace = TRUE)
                # Calculate bootstrap AUC
                return(sapply(preds, function(x) {
                        weighted.mean((x[boot_idx] - response[boot_idx])^2, weights[boot_idx])
                }))
        }
        
        return(t(replicate(B, run_boot())))
}

# Calibration by decile
set.seed(1)
boot_Brier_nhanes = boot_Brier_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

Brier_CI_df = cbind(
        "est" = Brier_nhanes, 
        "lci" = apply(boot_Brier_nhanes, 2, quantile, 0.025), 
        "uci" = apply(boot_Brier_nhanes, 2, quantile, 0.975))
Brier_CI_df

### ICI ###
calc_ici_weighted = function(preds, response, weights) {
        sapply(preds, function(x) {
                ICI = weighted.mean(abs(x - response), weights)
        })
}

ici_nhanes = calc_ici_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

# Function to calculate bootstrap weighted Brier score for multiple scores
boot_ici_weighted = function(preds, response,  weights, B = 1000) {
        run_boot = function(){
                # Sample indices with replacement
                boot_idx = sample(length(response), replace = TRUE)
                # Calculate bootstrap AUC
                return(sapply(preds, function(x) {
                        weighted.mean(abs(x[boot_idx] - response[boot_idx]), weights[boot_idx])
                }))
        }
        
        return(t(replicate(B, run_boot())))
}

set.seed(1)
boot_ici_nhanes = boot_ici_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

ici_CI_df = cbind(
        "est" = ici_nhanes, 
        "lci" = apply(boot_ici_nhanes, 2, quantile, 0.025), 
        "uci" = apply(boot_ici_nhanes, 2, quantile, 0.975))
ici_CI_df

### e50 ###
calc_e50_weighted = function(preds, response, weights) {
        sapply(preds, function(x) {
                abs_diff = abs(x - response)
                e50 = wtd.quantile(abs_diff, weights, probs = 0.50, type = "quantile")
        })
}

e50_nhanes = calc_e50_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

# Function to calculate bootstrap weighted Brier score for multiple scores
boot_e50_weighted = function(preds, response,  weights, B = 1000) {
        run_boot = function(){
                # Sample indices with replacement
                boot_idx = sample(length(response), replace = TRUE)
                # Calculate bootstrap AUC
                return(sapply(preds, function(x) {
                        wtd.quantile(abs(x[boot_idx] - response[boot_idx]), weights[boot_idx], probs = 0.50, type = "quantile")
                }))
        }
        
        return(t(replicate(B, run_boot())))
}

set.seed(1)
boot_e50_nhanes = boot_e50_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

e50_CI_df = cbind(
        "est" = e50_nhanes, 
        "lci" = apply(boot_e50_nhanes, 2, quantile, 0.025), 
        "uci" = apply(boot_e50_nhanes, 2, quantile, 0.975))
e50_CI_df


### e90 ###
calc_e90_weighted = function(preds, response, weights) {
        sapply(preds, function(x) {
                abs_diff = abs(x - response)
                e90 = wtd.quantile(abs_diff, weights, probs = 0.90, type = "quantile")
        })
}

e90_nhanes = calc_e90_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

# Function to calculate bootstrap weighted Brier score for multiple scores
boot_e90_weighted = function(preds, response,  weights, B = 1000) {
        run_boot = function(){
                # Sample indices with replacement
                boot_idx = sample(length(response), replace = TRUE)
                # Calculate bootstrap AUC
                return(sapply(preds, function(x) {
                        wtd.quantile(abs(x[boot_idx] - response[boot_idx]), weights[boot_idx], probs = 0.90, type = "quantile")
                }))
        }
        
        return(t(replicate(B, run_boot())))
}


set.seed(1)
boot_e90_nhanes = boot_e90_weighted(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP)

e90_CI_df = cbind(
        "est" = e90_nhanes, 
        "lci" = apply(boot_e90_nhanes, 2, quantile, 0.025), 
        "uci" = apply(boot_e90_nhanes, 2, quantile, 0.975))
e90_CI_df


# Function to calculate Cox intercept and slope for multiple scores, 
# while accounting for sampling weights
svy_cox_intercept_slope = function(preds, response, design) {
        sapply(preds, function(x) {
                my_df = data.frame(o = response, e = x)
                my_df$logite = log(my_df$e / (1 - my_df$e))
                design = svydesign(id = ~ nhanes$SDMVPSU,
                                   strata = ~ nhanes$SDMVSTRA, 
                                   weights = ~nhanes$WTSAFPRP, nest=TRUE, 
                                   data = my_df)
                
                fit = svyglm(o ~ I(logite), design = design, family = quasibinomial)
                est = coef(fit)
                ci = confint(fit)
                
                return(c(Intercept = paste0(round(est[1], 3), " (", round(ci[1,1], 3), 
                                            ", ", round(ci[1,2], 3), ")"), 
                         Slope = paste0(round(est[2], 3), " (", round(ci[2,1], 3), 
                                        ", ", round(ci[2,2], 3), ")")))
        })
}

# NHANES
svy_cox_intercept_slope(nhanes[, c("p_ccs", "eSL")], nhanes$fibrosis_f4_ft, nhanes$WTSAFPRP) %>% t()

# - Slope = direction of miscalibration, where 1 is perfect, greater than 1 means underestimating high risk/overestimating low risk, and lower than 1 means overestimating high risk/underestimating low risk 
# - Intercept = overall miscalibration, where 0 is good calibration, greater than 0 means average understimation, and lower than 0 means average overestimation. 
