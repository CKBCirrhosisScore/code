#####===Environment===#####
pacman::p_load(tidyverse,data.table,mice,ggpubr,MASS,rms,gridExtra,generalhoslem,pROC,knitr,dplyr,metamisc)

rm(list = ls())
gc()
setwd("")
load("Data/Model development.RData") 

#### =========== Penalized logistic regression ================ ####
f_full <- formula(fibrosis_f4 ~ rcs(age_3r, 3) + rcs(bmi_3r, 3) + rcs(whr_3r, 3) + 
                          rcs(ast_3r, 3) + rcs(ggt_3r, 3) + rcs(alt_3r, 3)+rcs(hdl_3r, 3) + rcs(tc_3r, 3) + rcs(tg_3r, 3)+ 
                          rcs(rpg_3r, 3) + rcs(sbp_3r, 3) + rcs(dbp_3r, 3) + is_female + region_is_urban + has_diabetes_3r + has_htn_3r +
                          is_female:(rcs(age_3r, 3) + rcs(bmi_3r, 3) + rcs(whr_3r, 3) + 
                                     rcs(ast_3r, 3) + rcs(ggt_3r, 3) +  rcs(alt_3r, 3)+
                                     rcs(hdl_3r, 3) + rcs(tc_3r, 3) + rcs(tg_3r, 3)+ rcs(rpg_3r, 3) + 
                                     rcs(sbp_3r, 3) + rcs(dbp_3r, 3) + is_female + region_is_urban + 
                                     has_diabetes_3r + has_htn_3r))

# Fit the maximal model (logistic regression)
m_full <- lrm(f_full, train_data)
a<- plot(anova(m_full), margin = c("d.f.", "chisq"))
plot(anova(m_full), margin = c("d.f.", "chisq"))

###*** Keep non-linear/interaction forms of variables according to chisq values ***###
###* Here non-linear forms of 1 variables (ggt_3r) according to chisq values and BMI*sex interaction as an example
f <- formula(fibrosis_f4 ~ age_3r  + whr_3r +
                     ast_3r + rcs(ggt_3r, 3) + alt_3r+ hdl_3r + tc_3r + tg_3r +
                     rpg_3r + sbp_3r + dbp_3r  + region_is_urban + has_diabetes_3r + has_htn_3r+
                     is_female*rcs(bmi_3r,3)) #+ rcs(bmi_3r,3)+ is_female
m <- lrm(f, train_data, x = TRUE, y = TRUE)
# Inspect fitted model
print(m, coefs = FALSE)

# We want the AIC to be lower for the reduced model than for the maximal model
# This indicates that the degrees of freedom we cut didn't lose us too much 
# predictive power
c(AIC_full = AIC(m_full), AIC_reduced = AIC(m))


# First find the optimal shrinkage parameter with a grid search
pen <- pentrace(m, list(simple = seq(0, 10, by = 1), 
                        nonlinear = seq(0, 10, by = 1),
                        interaction = seq(0, 10, by = 1)),which='aic')
opt.pen <- pen$penalty 

# One more iteration
pen <- pentrace(m, list(simple = do.call(seq, as.list(c(opt.pen$simple + c(0, 0.5), 0.1))), 
                     nonlinear = do.call(seq, as.list(c(opt.pen$nonlinear + c(0, 0.5), 0.1))),
                     interaction = do.call(seq, as.list(c(opt.pen$interaction + c(-0.5, 0.5), 0.1)))),,which='aic')

opt.pen <- pen$penalty

# Fit the model with penalized regression
m_pen <- update(m, penalty = list(simple = opt.pen$simple, 
                                  nonlinear = opt.pen$nonlinear,
                                  interaction = opt.pen$interaction))
print(m_pen, coefs = F)
AIC(m_pen)

# Predict from both the penalized and unpenalized model
train_data$lp0 <- predict(m, train_data, type = "lp")
train_data$p <- predict(m_pen, train_data, type = "fitted")
train_data$lp <- predict(m_pen, train_data, type = "lp")


# Find a simplified approximate model through ordinary linear regression on the 
# predicted values
(ols_f <- formula(paste("lp~", as.character(f))[3])) 
# lp ~ age_3r + whr_3r + ast_3r + rcs(ggt_3r, 3) + alt_3r + hdl_3r + 
#     tc_3r + tg_3r + rpg_3r + sbp_3r + dbp_3r + region_is_urban + 
#     has_diabetes_3r + has_htn_3r + is_female * rcs(bmi_3r, 3)

apx <- ols(ols_f, data = train_data, sigma = 1, x = TRUE)
fastbw(apx, aics = 1e10)
# 8 variables (ggt_3r, age_3r,  bmi_3r * is_female,  hdl_3r,has_diabetes_3r,  ast_3r ,tc_3r,is_female,bmi_3r ) are enough 
# to explain the predictions from the full model with R2>=0.95
tab2<- fastbw(apx, aics = 1e10)

###*** Keep only variables that make R2 at least 0.95 ***### ： 
#  8 variables (ggt_3r, age_3r,   hdl_3r,has_diabetes_3r,  ast_3r ,tc_3r,is_female,bmi_3r ) and 1 interaction (bmi_3r * is_female)
m_apx <- ols(fibrosis_f4 ~ rcs(ggt_3r, 3) +age_3r + tc_3r +is_female +  
                     ast_3r +  has_diabetes_3r +  hdl_3r +  rcs(bmi_3r,3) +
                     is_female*rcs(bmi_3r,3) , data = train_data, x = TRUE)

# Generate predicted values
train_data$lp2 <- predict(m_apx, train_data, type = "lp")
train_data$p2 <- plogis(train_data$lp2)

# Optionally force the approximated model into a logistic regression object
m_apx <- lrm(fibrosis_f4 ~ rcs(ggt_3r, 3) +age_3r + tc_3r +is_female +  ast_3r +  
has_diabetes_3r +  hdl_3r +  rcs(bmi_3r, 3) +is_female:rcs(bmi_3r, 3) , data = train_data, x = TRUE,
             penalty = list(simple = opt.pen$simple, 
                            nonlinear = opt.pen$nonlinear,
                            interaction = opt.pen$interaction))

# Replace the estimates with the appropriate ones based on the penalized model 
X <- cbind(Intercept = 1, m_pen$x) # full model design
Z <- cbind(Intercept = 1, m_apx$x) # approx. model design 
W <- solve(t(Z) %*% Z, t(Z)) %*% X   # contrast matrix
V <- vcov(m_pen) # var(full model)
V_apx <- W %*% V %*% t(W)

###adjusted, to fit the class type of original model###
coefficients1 <- W %*% m_pen$coefficients  #using contrast matrix of full model and final model to adjust the coefficients for final model
for(i in 1:length(m_apx$coefficients)){
        m_apx$coefficients[i] <- coefficients1[i]
}
m_apx$var <- V_apx

# Genereate predictions from the approximated model as a logistic model
train_data$lp2 <- predict(m_apx, train_data, type = "lp")
train_data$p2 <- predict(m_apx, train_data, type = "fitted")
auc(train_data$fibrosis_f4, train_data$p2)

m_apx
Function(m_apx) ###score
AIC(m_apx)

#### =========== eSL ================ ####
X <- train_data[, c("age_3r","whr_3r" ,"ast_3r" ,"ggt_3r", "alt_3r", "hdl_3r" , "tc_3r" , "tg_3r" ,
                    "rpg_3r" , "sbp_3r" , "dbp_3r"  , "region_is_urban" , "has_diabetes_3r" , "has_htn_3r",
                    "is_female","bmi_3r")]
Y4 <- train_data$fibrosis_f4

tune = list(ntrees = 500, max_depth = 1, shrinkage = 0.1)
XGB_learners = create.Learner("SL.xgboost", tune = tune, detailed_names = TRUE, name_prefix = "xgb")
SL_library_8 = c("SL.bayesglm", # Bayesian generalized linear model
                 "SL.earth", # Multivariate adaptive regression splines
                 "SL.gam", # Generalized additive model
                 "SL.gbm", # Generalized boosted model
                 "SL.glmnet", # Regularized generalized linear model
                 "SL.polymars", # Multivariate adaptive polynomial spline regression
                 "SL.randomForest", # Random forest
                 XGB_learners$names) # Support vector machine
Args <- base::commandArgs()
m <- as.numeric(Args[3]);m

set.seed(241227)
SL_8_F4 = SuperLearner(Y = Y4, X = X, family = binomial(), 
                       SL.library = SL_library_8, method = "method.AUC",
                       # stratified cross-validation scheme, 
                       # so that the outcome prevalence in training and validation sets 
                       # could be like the entire analytic dataset’s outcome prevalence
                       cvControl = list(V = 10, stratifyCV = TRUE))
  