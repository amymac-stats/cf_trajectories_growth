
source("./Data/Packages.R"); require(labelled)

base_uk  <- readRDS("./Data/base_uk.rds")
base_can  <- readRDS("./Data/base_can.rds")

all5 <- readRDS( file = "./Data/all5.RDS" )
cfc5 <- readRDS( file = "./Data/cfc5.RDS" )  

#----------------------#
# MISSING DATA DIAGRAM #
#----------------------#  


# INCLUDING: Tables at the end.


#-------#
# A. UK #
#-------#  


length(unique(all5$patient_id)) ## 3742 - this has had those with NO weights removed
nrow(base_uk) # 3758

process = data.frame( step = character(0), n = numeric(0) )
process[1, ] = c("Start", length(unique(all5$patient_id)) )


# STEPS

# 1. Number with a minimum age over 0.8 and maximum age over 5.5 

sum(base_uk$minAge<5.5) # 3758
sum(base_uk$maxAge>5.5) # 2373
table(base_uk$minAge<5.5 & base_uk$maxAge>5.5) # 2373

base3 <- subset(base_uk, minAge<5.5 & maxAge>5.5)
dim(base3) # 2373
table(base3$weightCount)  

process[2, ] = c("Within age range", nrow(base3) )


# 1. Number with at least one weight over the age of 0.8

sub1 = droplevels( subset(all5, age>=0.8 & maxAge>5.5 & !is.na(sds_weight_uk90) ) )
length(unique(sub1$patient_id)) # 2359

process[3, ] = c("At least one weight over 08", length(unique(sub1$patient_id)) )

# 2. At least one LF measure

base1 = subset(base_uk, patient_id %in% unique(sub1$patient_id))
base1 = subset(base1, !is.na(first_fev) & age_1st_fev<=7.5 &  !is.na(first_fvc) )
nrow(base1)

process[4, ] = c("At least one LF over 55 under 75", nrow(base1) )

# 3. No missing covariates 

sub2 <- subset(sub1, age>=0.8 &  !is.na(sds_weight_uk90) & !is.na(age) & !is.na(year_birth) &
                  !is.na(sex) & !is.na(dmg_IMD_zscore)   &
                  !is.na(pi) &   !is.na(F508_class) & 
                  !is.na(diagnosis)   &
                  !is.na(ageatdiagnosis)   &
                  !is.na(psa1) &  !is.na(staph1) &  !is.na(hflu1)  &
                  !is.na(SuppleFeed_any) )

base1$missing_covs <- ifelse( base1$patient_id %in% unique(sub2$patient_id), 
                                "Complete cases",
                                "Missing covariates" )

table(base1$missing_covs)

process[5, ] = c("No missing covariates", sum(base1$missing_covs=="Complete cases") )

# Summary

process$n = as.numeric(process$n)
process$diff = c(0, diff(process$n) )

# Marker for comparison table  

cc_ids_uk = base1 %>% filter(missing_covs=="Complete cases") %>% select(patient_id) %>% unlist()

base3$complete_case = base3$patient_id %in% cc_ids_uk

uk_all_cases = base3


#-----------#
# A. CANADA #
#-----------#  


length(unique(cfc5$patient_id)) ## 1562
nrow(base_can) # 1562

process_can = data.frame( step = character(0), n = numeric(0) )
process_can[1, ] = c("Start", length(unique(cfc5$patient_id)) )


# STEPS

# 1. Number with a minimum age over 0.8 and maximum age over 5.5 

sum(base_can$minAge<5.5) # 1562
sum(base_can$maxAge>5.5) # 979
table(base_can$minAge<5.5 & base_can$maxAge>5.5) # 979

base3 <- subset(base_can, minAge<5.5 & maxAge>5.5)
dim(base3) # 979
table(base3$weightCount)  

process_can[2, ] = c("Within age range", nrow(base3) )


# 1. Number with at least one weight over the age of 0.8

sub1 = droplevels( subset(cfc5, age>=0.8 & maxAge>5.5 & !is.na(sds_weight_uk90) ) )
length(unique(sub1$patient_id)) # 977

process_can[3, ] = c("At least one weight over 08", length(unique(sub1$patient_id)) )

# 2. At least one LF measure

base1 = subset(base_can, patient_id %in% unique(sub1$patient_id))
base1 = subset(base1, !is.na(first_fev) & age_1st_zfev<=7.5 &  !is.na(first_fvc) )
nrow(base1)

process_can[4, ] = c("At least one LF over 55 under 75", nrow(base1) )

# 3. No missing covariates 

sub2 <- subset(sub1, age>= 0.8 & !is.na(age) & !is.na(sex) & !is.na(pi) &  !is.na(F508_class) & 
                   !is.na(year_birth) & !is.na(diagnosis) & !is.na(ageatdiagnosis) & 
                   !is.na(psa1) & !is.na(hflu) & !is.na(staph1) )  
          

base1$missing_covs <- ifelse( base1$patient_id %in% unique(sub2$patient_id), 
                              "Complete cases",
                              "Missing covariates" )

table(base1$missing_covs)

process_can[5, ] = c("No missing covariates", sum(base1$missing_covs=="Complete cases") )

# Summary

process_can$n = as.numeric(process_can$n)
process_can$diff = c(0, diff(process_can$n) )


# Marker for comparison table  

cc_ids_can = base1 %>% filter(missing_covs=="Complete cases") %>% select(patient_id) %>% unlist()

base3$complete_case = base3$patient_id %in% cc_ids_can

can_all_cases = base3


# SAVE

write.csv(process, file = "./tables/missing_uk.csv")
write.csv(process_can, file = "./tables/missing_can.csv")



#-------------------#
# TABLES TO COMPARE #
#-------------------#  

# (i) Select variables

listVarsuk <- c( "sds_weight_uk90", "sds_bmi_uk90", "first_fev",
                 "sex", "year_birth", "pi","F508_class", "diagnosis", "ageatdiagnosis","weightCount", 
                 "dmg_IMD_zscore",
                 "Pa_sum", "Staph_sum", "Hflu_sum")

listVars_can <- c("sds_weight_uk90","sds_bmi_uk90", "first_fev", "sex", "year_birth", "pi",
                  "F508_class","diagnosis", "ageatdiagnosis", "weightCount",
                  "Pa_sum", "Staph_sum", "Hflu_sum") 

# (ii) Label then table variables

# UK

labelled::var_label(uk_all_cases$sex) <- "Sex"
labelled::var_label(uk_all_cases$year_birth) <- "Year of birth"
labelled::var_label(uk_all_cases$diagnosis) <- "Diagnosis"
labelled::var_label(uk_all_cases$ageatdiagnosis) <- "Age at diagnosis"
labelled::var_label(uk_all_cases$pi) <- "Pancreatic insufficient"
labelled::var_label(uk_all_cases$F508_class) <- "F508 genotype"
labelled::var_label(uk_all_cases$weightCount) <- "No. weight measurements"
labelled::var_label(uk_all_cases$dmg_IMD_zscore) <- "IMD z-score"

uk_all_cases$F508_class <- factor(uk_all_cases$F508_class, 
                             levels = levels(uk_all_cases$F508_class),
                             labels = c("Homozygous", "Heterozygous", "Other"))  


# CANADA

labelled::var_label(can_all_cases$sex) <- "Sex"
labelled::var_label(can_all_cases$year_birth) <- "Year of birth"
labelled::var_label(can_all_cases$diagnosis) <- "Diagnosis"
labelled::var_label(can_all_cases$ageatdiagnosis) <- "Age at diagnosis"
labelled::var_label(can_all_cases$pi) <- "Pancreatic insufficient"
labelled::var_label(can_all_cases$F508_class) <- "F508 genotype"
labelled::var_label(can_all_cases$weightCount) <- "No. weight measurements"


can_all_cases$diagnosis <- factor(can_all_cases$diagnosis, levels = levels(can_all_cases$diagnosis), labels = c("Asymptomatic", "MI", "Other", "Symptoms"))

can_all_cases$F508_class <- factor(can_all_cases$F508_class, levels = levels(can_all_cases$F508_class), labels = c("Homozygous", "Heterozygous", "Other"))


# (iii) Tables

tab_uk_missing =  KreateTableOne(x= uk_all_cases, 
                            vars = listVarsuk, 
                            strata = "complete_case",
                            factorVars=c("sex", "diagnosis","pi", "Genotype_class"),  
                            non_norm = c("year_birth", "dmg_IMD_zscore", 
                                         "ageatdiagnosis", "weightCount", "Pa_sum", "Staph_sum", "Hflu_sum"), digits = 1)


tab_can_missing = KreateTableOne(x = can_all_cases, 
                            vars = listVars_can, 
                            strata = "complete_case",
                            factorVars = c("sex", "diagnosis","pi", "Genotype_class"), 
                            test = TRUE, 
                            non_norm = c("year_birth", "ageatdiagnosis", "weightCount", "Pa_sum", "Staph_sum", "Hflu_sum"), digits = 1 )



# SAVE

write.csv(tab_uk_missing, file = "./tables/tab_uk_missing.csv")
write.csv(tab_can_missing, file = "./tables/tab_can_missing.csv")   
