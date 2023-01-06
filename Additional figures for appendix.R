
source("./Data/Packages.R"); require(labelled)


#-----------#
# LOAD DATA #
#-----------#  



base_uk  <- readRDS("./Data/base_uk.rds")
base_can  <- readRDS("./Data/base_can.rds")

all5 <- readRDS( file = "./Data/all5.RDS" )
cfc5 <- readRDS( file = "./Data/cfc5.RDS" )  



subuk <- subset(all5, age>=0.8 &  !is.na(sds_weight_uk90) & !is.na(age) & !is.na(year_birth) &
                  !is.na(sex) & !is.na(dmg_IMD_zscore)   &
                  !is.na(pi) &   !is.na(F508_class) & 
                  !is.na(diagnosis)   &
                  !is.na(ageatdiagnosis)   &
                  !is.na(psa1) &  !is.na(staph1) &  !is.na(hflu1)  &
                  !is.na(SuppleFeed_any) )  

subbase <- subset(base_uk, patient_id %in% unique(subuk$patient_id))

lf_df <- subset(subbase, !is.na(first_fev) & age_1st_fev<=7.5 &  !is.na(first_fvc) )


subcan <- subset(cfc5, age>= 0.8 & !is.na(age) & !is.na(sex) & !is.na(pi) &  !is.na(F508_class) & 
                   !is.na(year_birth) & !is.na(diagnosis) & !is.na(ageatdiagnosis) & 
                   !is.na(psa1) & !is.na(hflu) & !is.na(staph1) )

subbase_can <- subset(base_can, patient_id %in% unique(subcan$patient_id))

lf_dfc <- subset(subbase_can, !is.na(first_fev) & age_1st_zfev<=7.5 )


#----#
# Q1 #
#----#

# How many had only one weight/BMI z-score?

# UK
table(subbase$weightCount)
table(subbase$weightCount)/nrow(subbase)
# 293 (8.7%) of the UK had a single count

aggregate(age ~ weightCount, data = subbase, FUN = mean)
aggregate(year_birth ~ weightCount, data = subbase, FUN = median)


# Can
table(subbase_can$weightCount)
table(subbase_can$weightCount)/nrow(subbase_can)
# 47 (3.5%) in the Canadian data


#----#
# Q2 #
#----#  


# 3. While conceptualizing the initial weight z-score is straightforward, 
# it may be difficult to translate the difference in slopes into a clinical context. 
# The differences in slope per year seem quite small.

# Make it clear what a change of: 2.72 compared with a child with a slope of -0.08 would look like in real life.

# Idea: compare trajectories for a male/female who starts at the median weight z-score and has upper and lower quartile for growth

# starting with old code for plots
# aim is to translate into weight in kg, and also link to LF
# start by doing them all 

# Facet grid sex ~ measure (z-score, centile, kg)
# have the stable in solid lines
# and the decrease/increase in dashed lines
# label clearly

# Add to an appendix

#-----------------------------------------#
# LOAD DATA - with fitted values included #
#-----------------------------------------#


lf_df_plus = readRDS(file = "./Data/lf_df_plus.RDS")
subuk <- readRDS( file = "./Data/subuk.rds")


lf_dfc_plus = readRDS(file = "./Data/lf_dfc_plus.RDS")
subcan = readRDS( file = "./Data/subcan.rds")



#---------------------#
# INT + SLP QUANTILES #
#---------------------#

int_quant <- quantile( lf_df_plus$weight_intercept, c(0.25, 0.25,  0.50, 0.75, 0.75))
slp_quant <- quantile( lf_df_plus$weight_slope, c(0.75, 0.5, 0.50, 0.5, 0.25))/10 # this gives the value of the slope which matches the weight regression

eg_s <- data.frame( ints = int_quant,
                    slps = slp_quant,
                    int_label = c("low", "low", "mid", "high", "high"),
                    slp_label = c("high", "mid", "mid", "mid", "low"))


#-------------------------------#
# FIT AVERAGE WEIGHT TRAJECTORY #
#-------------------------------#

# for each combination of intercept and slope trajectory


fitted_mod <- readRDS( file = "./Objects/basic_mod.RDS" )

newdata = data.frame(age = seq(0,4.5, 0.1))
newdata$fits = predict(fitted_mod, newdata, re.form = NA)

newdata$low_high <- eg_s$slps[1]*newdata$age + newdata$fits + eg_s$ints[1]
newdata$low_mid  <- eg_s$slps[2]*newdata$age + newdata$fits + eg_s$ints[2]
newdata$mid_mid <- eg_s$slps[3]*newdata$age + newdata$fits + eg_s$ints[3]
newdata$high_mid <- eg_s$slps[4]*newdata$age + newdata$fits + eg_s$ints[4]
newdata$high_low <- eg_s$slps[5]*newdata$age + newdata$fits + eg_s$ints[5]


newdata <- gather(newdata, key = "lab", value = "fit", 3:7)

# Export to get weight in kg and weight centiles
# This is done in Excel

# newdata_m = newdata %>%
#   mutate(sex = 1) %>%
#   mutate(age = age + 1)
# 
# write.csv(newdata_m, "./Objects/weight_z_scores_m.csv", row.names = FALSE)


# newdata_f = newdata %>%
#   mutate(sex = 2) %>%
#   mutate(age = age + 1)
# 
# write.csv(newdata_f, "./Objects/weight_z_scores_f.csv", row.names = FALSE)

# Load to plot
ndf_m = read.csv(file = "./Objects/weight_z_scores_m_plus.csv") %>% mutate(sex = "Male")
ndf_f = read.csv(file = "./Objects/weight_z_scores_f_plus.csv") %>% mutate(sex = "Female")

ndf = rbind(ndf_m, ndf_f)

# Add variable for plotting

ndf$stable = ifelse( ndf$lab %in% c("high_mid",  "low_mid", "mid_mid"), "Stable", "Increase/decline" )



# Plot 


bigplot = ndf %>%
  select(age, sex, lab, stable, fit, centile_wt_uk, fitted_weight_kg) %>%
  rename(z_score = fit, centile = centile_wt_uk, kg = fitted_weight_kg) %>%
  pivot_longer(z_score:kg, names_to = "measure") %>%
  mutate( measure = fct_relevel(measure, "z_score", "centile", "kg")) %>%
  mutate( measure = fct_recode(measure, "z-score" = "z_score", "Centile" = "centile", "Kilograms" = "kg")) %>%
  mutate(Trajectory = fct_recode(lab, "High decline" = "high_low",
                          "High stable" = "high_mid",
                          "Low increase" = "low_high",
                          "Low stable" = "low_mid",
                          "Median stable" = "mid_mid")) %>%
  mutate(Trajectory = fct_relevel(Trajectory, "High stable", "High decline", "Median stable", "Low increase", "Low stable")) %>%
  ggplot(aes(age, value, colour = Trajectory)) +
  geom_line(aes(linetype = stable)) +
  scale_linetype_manual(values = c(2,1), name = NULL) +
  theme_bw(base_size = 12) +
  scale_colour_grey() +
  guides(linetype = "none") +
  labs(y = NULL, x = "Age, years") +
  theme(legend.position = "top") +
  facet_wrap(sex~measure, scales = "free_y" )  +
  guides(colour=guide_legend(nrow=2),byrow=TRUE) + # ,byrow=TRUE
  guides(linetype=guide_legend(nrow=2),byrow=TRUE)

ggsave(bigplot, file = "./Figures/For supplement/bigplot.png")
ggsave(bigplot, file = "./Figures/For supplement/bigplot.TIFF", width = 190, units = "mm", dpi = 1000)


females = ndf %>%
  filter(age %in% c(1,5)) %>%
  select(-fits) %>%
  rename(z_score = fit, centile = centile_wt_uk, kg = fitted_weight_kg) %>%
  arrange(age, lab, sex) %>%
  filter(sex=="Female")


males = ndf %>%
  filter(age %in% c(1,5)) %>%
  select(-fits) %>%
  rename(z_score = fit, centile = centile_wt_uk, kg = fitted_weight_kg) %>%
  arrange(age, lab, sex) %>%
  filter(sex=="Male")

# getting fitted lung functions

basic_mod_wt_fev = readRDS("./Objects/basic_mod_wt_fev.RDS")

df1 = eg_s %>% mutate(age_1st_fev = 6) %>% rename( slope = slps, intercept = ints) %>% mutate(slope = slope*10)
df1$pred = predict(basic_mod_wt_fev, newdata = df1, re.form = NA)
df1$lab = paste(df1$int_label, df1$slp_label, sep = "_")

to_merge = select(df1, lab, pred)

# merge in predicted FEV1%s

females = merge(females, to_merge, by = "lab") %>%
  mutate(lab = fct_recode(lab, "High decline" = "high_low",
                          "High stable" = "high_mid",
                          "Low increase" = "low_high",
                          "Low stable" = "low_mid",
                          "Median stable" = "mid_mid")) %>%
  mutate(lab = fct_relevel(lab, "High stable", "High decline", "Median stable", "Low increase", "Low stable")) %>%
  arrange(lab, age) %>%
  select(-stable) %>%
  rename("Predicted FEV1" = "pred")



males = merge(males, to_merge, by = "lab") %>%
  mutate(lab = fct_recode(lab, "High decline" = "high_low",
                          "High stable" = "high_mid",
                          "Low increase" = "low_high",
                          "Low stable" = "low_mid",
                          "Median stable" = "mid_mid")) %>%
  mutate(lab = fct_relevel(lab, "High stable", "High decline", "Median stable", "Low increase", "Low stable")) %>%
  arrange(lab, age) %>%
  select(-stable) %>%
  rename("Predicted FEV1" = "pred")


write.csv(females, file = "./tables/Tables for worked example/females.csv", row.names = FALSE)
write.csv(males, file = "./tables/Tables for worked example/males.csv", row.names = FALSE)


ndf %>%
  select(age, sex, lab, fit, SDS_Weight_UK_WHO) %>%
  rename(uk = fit, who = SDS_Weight_UK_WHO) %>%
  pivot_longer(uk:who, names_to = "ref", values_to = "z_score") %>%
  ggplot(aes(age, z_score, colour = lab)) +
  geom_line() +
  facet_wrap(ref~sex)






