

#------------#
# 2. TABLES  #
#------------#  

source("./Data/Packages.R"); require(labelled)

base_uk  <- readRDS("./Data/base_uk.rds")
base_can  <- readRDS("./Data/base_can.rds")

all5 <- readRDS( file = "./Data/all5.RDS" )
cfc5 <- readRDS( file = "./Data/cfc5.RDS" )

# EDIT 2020-12-16: table one for the main paper only includes those with a measure of LF. This does match the description in methods. 

#-----------------#
# 1. UK DATA PREP #
#-----------------#


subuk <- subset(all5, age>=0.8 &  !is.na(sds_weight_uk90) & !is.na(age) & !is.na(year_birth) &
                  !is.na(sex) & !is.na(dmg_IMD_zscore)   &
                  !is.na(pi) &   !is.na(F508_class) & 
                  !is.na(diagnosis)   &
                  !is.na(ageatdiagnosis)   &
                  !is.na(psa1) &  !is.na(staph1) &  !is.na(hflu1)  &
                  !is.na(SuppleFeed_any) )  

subbase <- subset(base_uk, patient_id %in% unique(subuk$patient_id))

### How many have lung function

lf_df <- subset(subbase, !is.na(first_fev) & age_1st_fev<=7.5 &  !is.na(first_fvc) )

base_uk$has_lf <- base_uk$patient_id %in% unique( lf_df$patient_id )

base_uk$missing_covs <- ifelse( base_uk$patient_id %in% unique(subuk$patient_id), 
                                "Complete cases",
                                "Missing covariates" )

listVarsuk <- c( "sds_weight_uk90", "sds_bmi_uk90", "first_fev", "first_zfev",
                 "sex", "year_birth", "pi","F508_class", "diagnosis", "ageatdiagnosis","weightCount", 
                 "dmg_IMD_zscore",
                 "Pa_sum", "Staph_sum", "Hflu_sum")



base_uk$F508_class <- factor(base_uk$F508_class, 
                             levels = levels(base_uk$F508_class),
                             labels = c("Homozygous", "Heterozygous", "Other"))  

base_uk$has_1_wt_fev <- as.vector(base_uk$has_1_wt_fev)


# label then table

labelled::var_label(base_uk$sex) <- "Sex"
labelled::var_label(base_uk$year_birth) <- "Year of birth"
labelled::var_label(base_uk$diagnosis) <- "Diagnosis"
labelled::var_label(base_uk$ageatdiagnosis) <- "Age at diagnosis"
labelled::var_label(base_uk$pi) <- "Pancreatic insufficient"
labelled::var_label(base_uk$F508_class) <- "F508 genotype"
labelled::var_label(base_uk$weightCount) <- "No. weight measurements"
labelled::var_label(base_uk$dmg_IMD_zscore) <- "IMD z-score"


#---------------------#
# 2. CANADA DATA PREP #
#---------------------#  


subcan <- subset(cfc5, age>= 0.8 & !is.na(age) & !is.na(sex) & !is.na(pi) &  !is.na(F508_class) & 
                   !is.na(year_birth) & !is.na(diagnosis) & !is.na(ageatdiagnosis) & 
                   !is.na(psa1) & !is.na(hflu) & !is.na(staph1) )



missing_vars <- subset(cfc5, age>= 0.8 | is.na(age) | is.na(sex) | is.na(pi) |  is.na(F508_class) | 
                         is.na(year_birth) | is.na(diagnosis) | is.na(ageatdiagnosis) | 
                         is.na(psa1) | is.na(hflu) | is.na(staph1) )



#---------#
# TABLE 1 #
#---------#

subcan <- subset(cfc5, age>= 0.8 & !is.na(age) & !is.na(sex) & !is.na(pi) &  !is.na(F508_class) & 
                   !is.na(year_birth) & !is.na(diagnosis) & !is.na(ageatdiagnosis) & 
                   !is.na(psa1) & !is.na(hflu) & !is.na(staph1) )

subbase_can <- subset(base_can, patient_id %in% unique(subcan$patient_id))

lf_dfc <- subset(subbase_can, !is.na(first_fev) & age_1st_zfev<=7.5 )

base_can$missing_covs <- ifelse( base_can$patient_id %in% unique(subcan$patient_id), 
                                 "Complete cases",
                                 "Missing covariates" )

base_can$has_lf <- base_can$patient_id %in% unique(lf_dfc$patient_id)

listVars_can <- c("sds_weight_uk90","sds_bmi_uk90", "first_fev", "first_zfev", "sex", "year_birth", "pi",
              "F508_class","diagnosis", "ageatdiagnosis", "weightCount",
              "Pa_sum", "Staph_sum", "Hflu_sum") 



# psa1 + staph1 + hflu1 +

# label then table

labelled::var_label(base_can$sex) <- "Sex"
labelled::var_label(base_can$year_birth) <- "Year of birth"
labelled::var_label(base_can$diagnosis) <- "Diagnosis"
labelled::var_label(base_can$ageatdiagnosis) <- "Age at diagnosis"
labelled::var_label(base_can$pi) <- "Pancreatic insufficient"
labelled::var_label(base_can$F508_class) <- "F508 genotype"
labelled::var_label(base_can$weightCount) <- "No. weight measurements"


base_can$diagnosis <- factor(base_can$diagnosis, levels = levels(base_can$diagnosis), labels = c("Asymptomatic", "MI", "Other", "Symptoms"))

base_can$F508_class <- factor(base_can$F508_class, levels = levels(base_can$F508_class), labels = c("Homozygous", "Heterozygous", "Other"))


#----------#
# 3. TABLE #
#----------#  

# For the main paper: ONLY THOSE with LF included


tab_uk_lf =  KreateTableOne(x= subset(base_uk, has_lf), 
                            vars = listVarsuk, 
                            factorVars=c("sex", "diagnosis","pi", "Genotype_class"),  
                            non_norm = c("year_birth", "dmg_IMD_zscore", 
                                         "ageatdiagnosis", "weightCount", "Pa_sum", "Staph_sum", "Hflu_sum"), digits = 1)


tab_can_lf = KreateTableOne(x = subset(base_can, has_lf), 
                         vars = listVars_can, 
                         factorVars = c("sex", "diagnosis","pi", "Genotype_class"), 
                         test = TRUE, 
                         non_norm = c("year_birth", "ageatdiagnosis", "weightCount", "Pa_sum", "Staph_sum", "Hflu_sum"), digits = 1 )




# SAVE

write.csv(tab_can_lf, file = "./tables/tab_can_lf.csv")
write.csv(tab_uk_lf, file = "./tables/tab_uk_lf.csv")   



