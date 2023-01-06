
#-------------#
# 5 - FIGURES #
#-------------#  

#----#
# UK #
#----#  

source("./Data/Packages.R"); require(labelled); require(ggrepel) 

base_uk  <- readRDS("./Data/base_uk.rds")
all5 <- readRDS( file = "./Data/all5.RDS" )


mod2 <- readRDS( file = "./Objects/mod1.RDS")
basic_mod <- readRDS(file = "./Objects/basicModuk.RDS")


subuk <- subset(all5, age>=0.8 &  !is.na(sds_weight_uk90) & !is.na(age) & !is.na(year_birth) &
                  !is.na(sex) & !is.na(dmg_IMD_zscore)   &
                  !is.na(pi) &   !is.na(F508_class) & 
                  !is.na(diagnosis)   &
                  !is.na(ageatdiagnosis)   &
                  !is.na(psa1) &  !is.na(staph1) &  !is.na(hflu1)  &
                  !is.na(SuppleFeed_any) )

subuk$age <- subuk$age - 1


lf_df <- subset(subbase, !is.na(first_fev) & age_1st_fev<=7.5 &  !is.na(first_fvc) )
 # added this Sept - should be fitting the weight model on only those with LF

lf_df_plus = readRDS(file = "./Data/lf_df_plus.RDS")
lf_dfc_plus = readRDS(file = "./Data/lf_dfc_plus.RDS")

subuk <- subset(subuk, patient_id %in% lf_df_plus$patient_id)
subcan <- subset(subcan, patient_id %in% lf_dfc_plus$patient_id)

#------------------------------------#
# SCATTERPLOT: LF by intercept/slope #
#------------------------------------#  



weight_intercept_fev = ggplot(lf_df_plus, aes(weight_intercept, first_fev)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") +
  geom_smooth(data = lf_dfc_plus, linetype = 3)

weight_slope_fev = ggplot(lf_df_plus, aes(weight_slope, first_fev)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 

bmi_intercept_fev = ggplot(lf_df_plus, aes(bmi_intercept, first_fev)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 

bmi_slope_fev = ggplot(lf_df_plus, aes(bmi_slope, first_fev)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 


grid.arrange(weight_intercept_fev, weight_slope_fev, bmi_intercept_fev, bmi_slope_fev, ncol = 2)




# FVC   



weight_intercept_fvc = ggplot(lf_df_plus, aes(weight_intercept, first_fvc)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 

weight_slope_fvc = ggplot(lf_df_plus, aes(weight_slope, first_fvc)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 

bmi_intercept_fvc = ggplot(lf_df_plus, aes(bmi_intercept, first_fvc)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 

bmi_slope_fvc = ggplot(lf_df_plus, aes(bmi_slope, first_fvc)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "red") 


grid.arrange(weight_intercept_fev, weight_slope_fev, bmi_intercept_fev, bmi_slope_fev, ncol = 2)
