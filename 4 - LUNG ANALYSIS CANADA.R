
#------------------------------#
# LUNG FUNCTION MODELS: CANADA #
#------------------------------#  

# 10/9/21: Fit two models (i) With no covariates other than slope/intercept; (ii) with only baseline covariates. 
# 23/9/21: Model fitting for weight/BMI model - testing for non-linear terms for age


#------#
# DATA #
#------#   


source("./Data/Packages.R"); require(labelled)

lf_dfc = readRDS( file = "./Data/lf_dfc.rds")
subcan = readRDS( file = "./Data/subcan.rds")

subcan[subcan$patient_id=="197" & subcan$age>4, "sds_bmi_uk90"] <- NA # 4/10/21 - found this error (BMI of 2)


#--------------#
# 0. FUNCTIONS #
#--------------#


boot_fun_uk <- function( seed = 1, df1 = subcan, df2 = base_can, id = "patient",
                         wt_outcome = "sds_weight_uk90",
                         lf_formula = "first_fev ~ age_1st_fev + 
                         intercept + 
                         slope",
                         return_data = TRUE, 
                         return_mod = FALSE){
  
  if(length(df2$intercept)>0){ # the intercept will be there iff the rest are
    df2$int_raw <- NULL        # removing the present values of intercept etc.
    df2$intercept <- NULL
    df2$slp_raw <- NULL
    df2$slope <- NULL
  }
  
  
  if(!is.na(seed)){ 
    set.seed(seed)
    
    ids <- droplevels(unique(df2$patient)) # same from the baseline dataset
    samp <- sample(ids, length(ids), replace = TRUE) # THE SAMPLE  
    
    # BASE DATA
    
    reps <- data.frame( table(samp) )
    subreps <- subset( reps, samp %in% df2$patient )
    names( subreps ) <- c("patient_id", "reps")
    
    df2 <- merge( df2, subreps, by = "patient_id" ) # how many times should each ID be repeated
    df2 <- df2[ df2$rep!=0, ] # remove all those with zero reps
    df2 <- df2[ rep(1:dim(df2)[1], df2$rep), ] # rep each row number of bootstrap replications 
    
    dfb <-  df2      # BOOTSTAPPED BASE DATA
    
  } else {
    
    # to get the coefficient point estimates from the original data
    dfb <- df2
  }
  
  # Model fitting for weight/BMI model - testing for non-linear terms for age
  
  formulae = list( "~  age + (age|patient)" = as.formula( paste( wt_outcome, "~  age + (age|patient)") ),
                   "~  age + ns(age, 2) + (age|patient)" = as.formula( paste( wt_outcome, "~  age + ns(age, 2) + (age|patient)") ),
                   "~  age + ns(age, 3) + (age|patient)" = as.formula( paste( wt_outcome, "~  age + ns(age, 3) + (age|patient)") ) )
  
  chosen_formula <- formulae %>%
    map( ~ lmer(. , data = df1, control = lmerControl(optimizer ="Nelder_Mead") ) ) %>%
    map( ~ BIC(.)) %>%
    unlist() %>%
    abs() %>%
    which.min() %>%
    names()
  
  chosen_formula = as.formula( paste( wt_outcome, chosen_formula ) )
  
  print(paste0("Number of knots: ", which(formulae == chosen_formula)) )
  
  mod_basic = lmer(chosen_formula , data = df1, control = lmerControl(optimizer ="Nelder_Mead") )
  
  # estimate the ranefs
  ranfs <- ranef(mod_basic)[[id]]
  
  ranfs[[id]] <- rownames(ranfs)
  names(ranfs)[1] <- "intercept"
  names(ranfs)[2] <- "slope"
  
  dfb <- merge(dfb, ranfs, by = id) # merge with the re-sampled baseline data  
  
  # divide ranefs by 10 - ONLY THE SLOPE?
  # but also return the originals
  dfb$int_raw <- dfb$intercept
  dfb$slp_raw <- dfb$slope
  
  dfb$intercept <- dfb$int_raw
  dfb$slope <- dfb$slp_raw*10
  
  # edit year of birth
  dfb$year_birth <- dfb$year_birth - 2000
  
  # lung function linear regression
  fevmod3 <- lm( lf_formula,
                 data = dfb )
  
  if(!is.na(seed) | !return_data){ 
    coefs <- coef(fevmod3)
  }
  
  if(is.na(seed) & return_data ){
    coefs <- list( betas = coef(fevmod3), data = dfb)
  }
  if(return_mod == TRUE){
    coefs$model <- fevmod3 
  }
  
  return(coefs) 
  
}



#-------------#
# 1. FORMULAE #
#-------------#  

wt_formula <- "sds_weight_uk90"
bmi_formula <- "sds_bmi_uk90"

fev_formula_basic = "first_fev ~ age_1st_zfev + 
                         intercept + 
                         slope"

fvc_formula_basic = "first_fvc ~ age_1st_zfev + 
                         intercept + 
                         slope"


fev_formula_adj = "first_fev ~ age_1st_zfev + 
                         intercept + 
                         slope + sex + year_birth +
                         F508_class + pi + 
                         diagnosis +
                         ageatdiagnosis"

fvc_formula_adj = "first_fvc ~ age_1st_zfev + 
                         intercept + 
                         slope + sex + year_birth +
                         F508_class + pi + 
                         diagnosis +
                         ageatdiagnosis"



#--------------#
# 1.(b) LABELS #
#--------------#

labs_basic_wt = c("Intercept", "Age at lung function test",  
                  "Weight z-score age 1\nunits = 1 z-score", 
                  "Weight z-score slope\nunits = 0.1 z-scores")

labs_basic_bmi = c("Intercept", "Age at lung function test",  
                   "BMI z-score age 1\nunits = 1 z-score", 
                   "BMI z-score slope\nunits = 0.1 z-scores")


labs_adj_wt = c("Intercept", "Age at lung function test",  
                "Weight z-score age 1\nunits = 1 z-score", 
                "Weight z-score slope\nunits = 0.1 z-scores",
                "Sex = FEMALE", "Year of birth", 
                "F508 Heterozygous", "F508 Other", 
                "PI",
                "Diagnosis: MI",
                "Diagnosis: Other",
                "Diagnosis: Symptoms", 
                "Age at diagnosis")

labs_adj_bmi = c("Intercept", "Age at lung function test",  
                 "BMI z-score age 1\nunits = 1 z-score", 
                 "BMI z-score slope\nunits = 0.1 z-scores",
                 "Sex = FEMALE", "Year of birth", 
                 "F508 Heterozygous", "F508 Other", 
                 "PI",
                 "Diagnosis: MI",
                 "Diagnosis: Other",
                 "Diagnosis: Symptoms", 
                 "Age at diagnosis")


#----------------------#
# 2 UNADJUSTED MODELS  #
#----------------------#



#------------#
# (a) WEIGHT #
#------------#

### (i) FEV1%

# Model
basic_mod_wt_fev <- boot_fun_uk(seed = NA,
                                wt_outcome = wt_formula, 
                                lf_formula = fev_formula_basic,
                                df1 = subcan,
                                df2 = lf_dfc, 
                                return_data = FALSE, return_mod = TRUE)$model


# Boots
# booted_coefs_wt1 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                         wt_outcome = wt_formula,
#                                                                         lf_formula = fev_formula_basic,
#                                                                         df1 = subcan,
#                                                                      df2 = lf_dfc) ))
# 
# saveRDS(booted_coefs_wt1, file = "./Objects/booted_coefs_can_wt_1.RDS")

booted_coefs_wt1 <- readRDS( file = "./Objects/booted_coefs_can_wt_1.RDS" )

# CI
confints_wt1 <- data.frame( t(apply(booted_coefs_wt1, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_wt1) <- c("boot_lwr", "boot_upr")  


# Save results
coef_wt_fev1 = coef(basic_mod_wt_fev)
coefs_wt1 <- data.frame( coef = coef_wt_fev1, lwr = confints_wt1$boot_lwr, upr = confints_wt1$boot_upr )

coefs_wt1$label = labs_basic_wt
coefs_wt1$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_wt1$coef, coefs_wt1$lwr, coefs_wt1$upr)

write.csv(coefs_wt1, file = "./Objects/Model_results/can_confints_wt_fev.csv")





### (ii) FVC%  

# Model
basic_mod_wt_fvc <- boot_fun_uk(seed = NA,
                                wt_outcome = wt_formula, 
                                lf_formula = fvc_formula_basic,
                                df1 = subcan,
                                df2 = lf_dfc, 
                                return_data = FALSE, return_mod = TRUE)$model

# Boots
# booted_coefs_wt_fvc1 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                        wt_outcome = wt_formula,
#                                                                        lf_formula = fvc_formula_basic,
#                                                                        df1 = subcan,
#                                                                     df2 = lf_dfc) )) #
# 
# saveRDS(booted_coefs_wt_fvc1, file = "./Objects/booted_coefs_wt_fvc1.RDS")

booted_coefs_wt_fvc1 <- readRDS( file = "./Objects/booted_coefs_wt_fvc1.RDS" )

# CI
confints_wt_fvc1 <- data.frame( t(apply(booted_coefs_wt_fvc1, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_wt_fvc1) <- c("boot_lwr", "boot_upr") 


# Save results
coef_wt_fvc1 = coef(basic_mod_wt_fvc)
coefs_wt_fvc1 <- data.frame( coef = coef_wt_fvc1, lwr = confints_wt_fvc1$boot_lwr, upr = confints_wt_fvc1$boot_upr )

coefs_wt_fvc1$label = labs_basic_wt
coefs_wt_fvc1$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_wt_fvc1$coef, coefs_wt_fvc1$lwr, coefs_wt_fvc1$upr)

write.csv(coefs_wt_fvc1, file = "./Objects/Model_results/can_confints_wt_fvc.csv")


#---------#
# (b) BMI #
#---------#  

### (i) FEV1%


# Model
basic_mod_bmi_fev <- boot_fun_uk(seed = NA, # didn't converge
                                 wt_outcome = bmi_formula, 
                                 lf_formula = fev_formula_basic,
                                 df1 = subcan,
                                 df2 = lf_dfc, 
                                 return_data = FALSE, return_mod = TRUE)$model

# Boots
# booted_coefs_bmi1 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                         wt_outcome = bmi_formula,
#                                                                        lf_formula = fev_formula_basic,
#                                                                        df1 = subcan,
#                                                                     df2 = lf_dfc) )) #

saveRDS(booted_coefs_bmi1, file = "./Objects/booted_coefs_can_bmi_1.RDS")

booted_coefs_bmi1 <- readRDS( file = "./Objects/booted_coefs_can_bmi_1.RDS" )


# CI
confints_bmi1 <- data.frame( t(apply(booted_coefs_bmi1, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_bmi1) <- c("boot_lwr", "boot_upr")


# Save results
coef_bmi_fev1 = coef(basic_mod_bmi_fev)
coefs_bmi_fev1 <- data.frame( coef = coef_bmi_fev1, lwr = confints_bmi1$boot_lwr, upr = confints_bmi1$boot_upr )

coefs_bmi_fev1$label = labs_basic_bmi
coefs_bmi_fev1$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_bmi_fev1$coef, coefs_bmi_fev1$lwr, coefs_bmi_fev1$upr)

write.csv(coefs_bmi_fev1, file = "./Objects/Model_results/can_confints_bmi_fev.csv")





### (ii) FVC%  

# Model
basic_mod_bmi_fvc <- boot_fun_uk(seed = NA, # didnt converge
                                 wt_outcome = bmi_formula, 
                                 lf_formula = fvc_formula_basic,
                                 df1 = subcan,
                                 df2 = lf_dfc, 
                                 return_data = FALSE, return_mod = TRUE)$model

# Boots
# booted_coefs_bmi_fvc1 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                             wt_outcome = bmi_formula,
#                                                                            lf_formula = fvc_formula_basic,
#                                                                            df1 = subcan,
#                                                                            df2 = lf_dfc) )) #
# 
# saveRDS(booted_coefs_bmi_fvc1, file = "./Objects/booted_coefs_bmi_fvc1.RDS")

booted_coefs_bmi_fvc1 <- readRDS( file = "./Objects/booted_coefs_bmi_fvc1.RDS" )

# CI
confints_bmi_fvc1 <- data.frame( t(apply(booted_coefs_bmi_fvc1, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_bmi_fvc1) <- c("boot_lwr", "boot_upr")

# Save results
coef_bmi_fvc1 = coef(basic_mod_bmi_fvc)
coefs_bmi_fvc1 <- data.frame( coef = coef_bmi_fvc1, lwr = confints_bmi_fvc1$boot_lwr, upr = confints_bmi_fvc1$boot_upr )

coefs_bmi_fvc1$label = labs_basic_bmi
coefs_bmi_fvc1$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_bmi_fvc1$coef, coefs_bmi_fvc1$lwr, coefs_bmi_fvc1$upr)

write.csv(coefs_bmi_fvc1, file = "./Objects/Model_results/can_confints_bmi_fvc.csv")



## there's more... ##








#--------------------#
# 3 ADJUSTED MODELS  #
#--------------------#



#------------#
# (a) WEIGHT #
#------------#

### (i) FEV1%

# Model
adj_mod_wt_fev <- boot_fun_uk(seed = NA,
                              wt_outcome = wt_formula, 
                              lf_formula = fev_formula_adj,
                              df1 = subcan,
                              df2 = lf_dfc, 
                              return_data = FALSE, return_mod = TRUE)$model


# Boots
# booted_coefs_wt2 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                        wt_outcome = wt_formula,
#                                                                        lf_formula = fev_formula_adj,
#                                                                        df1 = subcan,
#                                                                        df2 = lf_dfc) ))
# 
# saveRDS(booted_coefs_wt2, file = "./Objects/booted_coefs_can_wt_2.RDS")

booted_coefs_wt2 <- readRDS( file = "./Objects/booted_coefs_can_wt_2.RDS" )

# CI
confints_wt2 <- data.frame( t(apply(booted_coefs_wt2, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_wt2) <- c("boot_lwr", "boot_upr")  


# Save results
coef_wt_fev2 = coef(adj_mod_wt_fev)
coefs_wt2 <- data.frame( coef = coef_wt_fev2, lwr = confints_wt2$boot_lwr, upr = confints_wt2$boot_upr )

coefs_wt2$label = labs_adj_wt
coefs_wt2$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_wt2$coef, coefs_wt2$lwr, coefs_wt2$upr)

write.csv(coefs_wt2, file = "./Objects/Model_results/can_confints_wt_fev2.csv")





### (ii) FVC%  

# Model
adj_mod_wt_fvc <- boot_fun_uk(seed = NA,
                              wt_outcome = wt_formula, 
                              lf_formula = fvc_formula_adj,
                              df1 = subcan,
                              df2 = lf_dfc, 
                              return_data = FALSE, return_mod = TRUE)$model

# Boots
# booted_coefs_wt_fvc2 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                            wt_outcome = wt_formula,
#                                                                            lf_formula = fvc_formula_adj,
#                                                                            df1 = subcan,
#                                                                            df2 = lf_dfc) )) #
# 
# saveRDS(booted_coefs_wt_fvc2, file = "./Objects/booted_coefs_wt_fvc2.RDS")

booted_coefs_wt_fvc2 <- readRDS( file = "./Objects/booted_coefs_wt_fvc2.RDS" )

# CI
confints_wt_fvc2 <- data.frame( t(apply(booted_coefs_wt_fvc2, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_wt_fvc2) <- c("boot_lwr", "boot_upr") 


# Save results
coef_wt_fvc2 = coef(adj_mod_wt_fvc)
coefs_wt_fvc2 <- data.frame( coef = coef_wt_fvc2, lwr = confints_wt_fvc2$boot_lwr, upr = confints_wt_fvc2$boot_upr )

coefs_wt_fvc2$label = labs_adj_wt
coefs_wt_fvc2$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_wt_fvc2$coef, coefs_wt_fvc2$lwr, coefs_wt_fvc2$upr)

write.csv(coefs_wt_fvc2, file = "./Objects/Model_results/can_confints_wt_fvc2.csv")


#---------#
# (b) BMI #
#---------#  

### (i) FEV1%


# Model
adj_mod_bmi_fev <- boot_fun_uk(seed = NA, # didnt converge
                               wt_outcome = bmi_formula, 
                               lf_formula = fev_formula_adj,
                               df1 = subcan,
                               df2 = lf_dfc, 
                               return_data = FALSE, return_mod = TRUE)$model

# Boots
# booted_coefs_bmi2 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                         wt_outcome = bmi_formula,
#                                                                         lf_formula = fev_formula_adj,
#                                                                         df1 = subcan,
#                                                                         df2 = lf_dfc) )) #
# 
# saveRDS(booted_coefs_bmi2, file = "./Objects/booted_coefs_can_bmi_2.RDS")

booted_coefs_bmi2 <- readRDS( file = "./Objects/booted_coefs_can_bmi_2.RDS" )


# CI
confints_bmi2 <- data.frame( t(apply(booted_coefs_bmi2, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_bmi2) <- c("boot_lwr", "boot_upr")


# Save results
coef_bmi_fev2 = coef(adj_mod_bmi_fev)
coefs_bmi_fev2 <- data.frame( coef = coef_bmi_fev2, lwr = confints_bmi2$boot_lwr, upr = confints_bmi2$boot_upr )

coefs_bmi_fev2$label = labs_adj_bmi
coefs_bmi_fev2$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_bmi_fev2$coef, coefs_bmi_fev2$lwr, coefs_bmi_fev2$upr)

write.csv(coefs_bmi_fev2, file = "./Objects/Model_results/can_confints_bmi_fev2.csv")




### (ii) FVC%  

# Model
adj_mod_bmi_fvc <- boot_fun_uk(seed = NA, # didn't converge
                               wt_outcome = bmi_formula, 
                               lf_formula = fvc_formula_adj,
                               df1 = subcan,
                               df2 = lf_dfc, 
                               return_data = FALSE, return_mod = TRUE)$model

# Boots
# booted_coefs_bmi_fvc2 <- data.frame(sapply( 1:1000, function(x) boot_fun_uk(seed = x,
#                                                                             wt_outcome = bmi_formula,
#                                                                             lf_formula = fvc_formula_adj,
#                                                                             df1 = subcan,
#                                                                             df2 = lf_dfc) )) #
# 
# saveRDS(booted_coefs_bmi_fvc2, file = "./Objects/booted_coefs_bmi_fvc2.RDS")

booted_coefs_bmi_fvc2 <- readRDS( file = "./Objects/booted_coefs_bmi_fvc2.RDS" )

# CI
confints_bmi_fvc2 <- data.frame( t(apply(booted_coefs_bmi_fvc2, 1, function(x) quantile(x, probs = c(0.025, 0.975))) ))
names(confints_bmi_fvc2) <- c("boot_lwr", "boot_upr")

# Save results
coef_bmi_fvc2 = coef(adj_mod_bmi_fvc)
coefs_bmi_fvc2 <- data.frame( coef = coef_bmi_fvc2, lwr = confints_bmi_fvc2$boot_lwr, upr = confints_bmi_fvc2$boot_upr )

coefs_bmi_fvc2$label = labs_adj_bmi
coefs_bmi_fvc2$final <- sprintf("%.2f [%.2f, %.2f]" , coefs_bmi_fvc2$coef, coefs_bmi_fvc2$lwr, coefs_bmi_fvc2$upr)

write.csv(coefs_bmi_fvc2, file = "./Objects/Model_results/can_confints_bmi_fvc2.csv")







#------------------------------------#
# GET PREDICTED WEIGHT/BMI FOR PLOTS #
#------------------------------------#  

# newdat = data.frame( age = seq(0,4.5,0.5) )
# wt_mod = lmer(wt_outcome , data = subuk )
# newdat$weight_preds = predict(wt_mod, newdata = newdat, re.form = NA )

# 4/10/21 - want to add predicted LF
# select the person with median, upper and lower quartile intercepts, take whatever slope they actually have,
# then find predicted LF



weight_summaries = boot_fun_uk(seed = NA,
                               wt_outcome = wt_formula, 
                               lf_formula = fev_formula_adj,
                               df1 = subcan,
                               df2 = lf_dfc, 
                               return_data = TRUE, return_mod = FALSE)$data %>%
  select(patient_id, intercept, slope) %>%
  rename(weight_intercept = intercept, weight_slope = slope)



bmi_summaries = boot_fun_uk(seed = NA,
                            wt_outcome = bmi_formula, 
                            lf_formula = fev_formula_adj,
                            df1 = subcan,
                            df2 = lf_dfc, 
                            return_data = TRUE, return_mod = FALSE)$data %>%
  select(patient_id, intercept, slope) %>%
  rename(bmi_intercept = intercept, bmi_slope = slope)




# merge all  


lf_dfc_plus = merge(lf_dfc, weight_summaries, by = "patient_id")
lf_dfc_plus = merge(lf_dfc_plus, bmi_summaries, by = "patient_id")


# predicted FEV and FVC  

# Divide people into thirds according to intercept and predict their lung function.

wt_lf_preds = lf_dfc_plus %>%
  select(-bmi_intercept, -bmi_slope) %>%
  rename( intercept = weight_intercept, slope = weight_slope ) %>%
  mutate( predicted_wt_fev = predict(basic_mod_wt_fev, newdata = .),
          predicted_wt_fvc = predict(basic_mod_wt_fvc, newdata = .)) %>%
  select(patient_id, predicted_wt_fev, predicted_wt_fvc)  


bmi_lf_preds = lf_dfc_plus %>%
  select(-weight_intercept, -weight_slope) %>%
  rename( intercept = bmi_intercept, slope = bmi_slope ) %>%
  mutate( predicted_bmi_fev = predict(basic_mod_bmi_fev, newdata = .),
          predicted_bmi_fvc = predict(basic_mod_bmi_fvc, newdata = .)) %>%
  select(patient_id, predicted_bmi_fev, predicted_bmi_fvc)


# merge again

lf_dfc_plus = merge(lf_dfc_plus, wt_lf_preds, by = "patient_id")
lf_dfc_plus = merge(lf_dfc_plus, bmi_lf_preds, by = "patient_id")


# save

saveRDS(lf_dfc_plus, file = "./Data/lf_dfc_plus.RDS")




#
# Extra bits 
#

df1 = data.frame( age_1st_zfev = 6, slope = quantile(lf_dfc_plus$weight_slope, c(0.25, 0.75)), intercept = 0)

df1$pred = predict(basic_mod_wt_fev, newdata = df1, re.form = NA)
df1$pred[2] - df1$pred[1]        


df2 = data.frame( age_1st_zfev = 6, intercept = quantile(lf_dfc_plus$weight_intercept, c(0.25, 0.75)), slope = 0)       

df2$pred = predict(basic_mod_wt_fev, newdata = df2, re.form = NA)
df2$pred[2] - df2$pred[1]  

saveRDS(df1, "./Objects/df1_for_scatter_can.RDS")
saveRDS(df2, "./Objects/df2_for_scatter_can.RDS")

 
