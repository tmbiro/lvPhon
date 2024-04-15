#command + option + T = run entire section
#command + shift + R = add new section
#command + option + O = fold all sections
#command + option + shift + O = unfold all sections
#command + enter = run line

# Set-up environment ------------------------------------------------------------------

# Clear workspace
rm(list = ls(all = TRUE))

library(broom)
library(R.utils)
library(lmerTest)
library(LMERConvenienceFunctions)
library(car)
library(lme4)
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)

# Set options 
options(max.print=100000)

# Define helper functions
## For mean centering
myCenter = function(x) x - mean(x, na.rm=T)

## For logging values
Logit = function(x) log(x/(1-x))

# Set up datasets & check distributions ---------------------------------------------------------

#setwd('/Users/tifanibiro/Library/CloudStorage/GoogleDrive-tifbiro@gmail.com/My Drive/Research/Aphasia/lvPhon/Scripting/Brendan/Rebased/lvPhon/R/Resources')
setwd('H:/My Drive/Research/Aphasia/lvPhon/Scripting/Brendan/Rebased/lvPhon/R/Resources')
data <- as.data.frame(read.csv(file='lvPhon-MRI-Acc-dat-v4.csv', header=TRUE, sep=','))
head(data)

# Make a copy of the data file for volume 
volDat <- subset(
  data, 
  measure == 'volume', 
  select = c(
    PID,
    Region_ID_Num,
    measure,
    metric,
    value,
    Phon_Acc_Group,
    Place_Acc_Group,
    Manner_Acc_Group,
    Voicing_Acc_Group,
    Phon_Acc,
    Manner_Acc,
    Place_Acc,
    Voicing_Acc,
    Lobe,
    Gyrus_Label,
    Gyrus,
    Hemisphere
    ))

head(volDat)

# Declare factors for volumetric dataset
volDat$Hemisphere <- ifelse(volDat$Hemisphere == 'L', 'Left', 'Right')
volDat$Hemisphere <- factor(volDat$Hemisphere, levels = c('Right','Left'))
volDat$PID <- factor(volDat$PID)
volDat$Region_ID_Num <- factor(volDat$Region_ID_Num)
volDat$measure <- factor(volDat$measure)
volDat$metric <- factor(volDat$metric)
volDat$Gyrus <- factor(volDat$Gyrus)
volDat$Gyrus_Label <- factor(volDat$Gyrus_Label)
volDat$Lobe <- factor(volDat$Lobe)
volDat$Phon_Acc <- scale(volDat$Phon_Acc)
volDat$Manner_Acc <- scale(volDat$Manner_Acc)
volDat$Place_Acc <- scale(volDat$Place_Acc)
volDat$Voicing_Acc <- scale(volDat$Voicing_Acc)
volDat$Phon_Acc_Group <- factor(volDat$Phon_Acc_Group, 
                               levels = c(
                                 '96-100%',
                                 '91-95%',
                                 '86-90%',
                                 '81-85%',
                                 '76-80%',
                                 '71-75%',
                                 '66-70%', 
                                 '61-65%'
                                   ))
volDat$Place_Acc_Group <- factor(volDat$Place_Acc_Group, 
                                levels = c(
                                  '96-100%',
                                  '91-95%',
                                  '86-90%',
                                  '81-85%',
                                  '76-80%',
                                  '71-75%'
                                ))
volDat$Voicing_Acc_Group <- factor(volDat$Voicing_Acc_Group, 
                                 levels = c(
                                   '96-100%',
                                   '91-95%',
                                   '86-90%'
                                 ))
volDat$Manner_Acc_Group <- factor(volDat$Manner_Acc_Group, 
                                levels = c(
                                  '96-100%',
                                  '91-95%',
                                  '86-90%',
                                  '81-85%',
                                  '76-80%',
                                  '71-75%'
                                ))
summary(volDat)

# Make a copy of the data file for thickness 
thicDat <- subset(
  data, 
  measure == 'thickness', 
  select = c(
    PID,
    Region_ID_Num,
    measure,
    metric,
    value,
    Phon_Acc_Group,
    Place_Acc_Group,
    Manner_Acc_Group,
    Voicing_Acc_Group,
    Phon_Acc,
    Manner_Acc,
    Place_Acc,
    Voicing_Acc,
    Lobe,
    Gyrus_Label,
    Gyrus,
    Hemisphere
  ))

# Take out only means
thicDat2 <- subset(thicDat, metric == 'mean')

head(thicDat2)

# Declare factors for thickness dataset
thicDat2$Hemisphere <- ifelse(thicDat2$Hemisphere == 'L', 'Left', 'Right')
thicDat2$Hemisphere <- factor(thicDat2$Hemisphere, levels = c('Right','Left'))
thicDat2$PID <- factor(thicDat2$PID)
thicDat2$Region_ID_Num <- factor(thicDat2$Region_ID_Num)
thicDat2$measure <- factor(thicDat2$measure)
thicDat2$metric <- factor(thicDat2$metric)
thicDat2$Gyrus <- factor(thicDat2$Gyrus)
thicDat2$Gyrus_Label <- factor(thicDat2$Gyrus_Label)
thicDat2$Lobe <- factor(thicDat2$Lobe)
thicDat2$Phon_Acc <- scale(thicDat2$Phon_Acc)
thicDat2$Manner_Acc <- scale(thicDat2$Manner_Acc)
thicDat2$Place_Acc <- scale(thicDat2$Place_Acc)
thicDat2$Voicing_Acc <- scale(thicDat2$Voicing_Acc)
thicDat2$Phon_Acc_Group <- factor(thicDat2$Phon_Acc_Group, 
                                levels = c(
                                  '96-100%',
                                  '91-95%',
                                  '86-90%',
                                  '81-85%',
                                  '76-80%',
                                  '71-75%',
                                  '66-70%', 
                                  '61-65%'
                                ))
thicDat2$Place_Acc_Group <- factor(thicDat2$Place_Acc_Group, 
                                 levels = c(
                                   '96-100%',
                                   '91-95%',
                                   '86-90%',
                                   '81-85%',
                                   '76-80%',
                                   '71-75%'
                                 ))
thicDat2$Voicing_Acc_Group <- factor(thicDat2$Voicing_Acc_Group, 
                                   levels = c(
                                     '96-100%',
                                     '91-95%',
                                     '86-90%'
                                   ))
thicDat2$Manner_Acc_Group <- factor(thicDat2$Manner_Acc_Group, 
                                  levels = c(
                                    '96-100%',
                                    '91-95%',
                                    '86-90%',
                                    '81-85%',
                                    '76-80%',
                                    '71-75%'
                                  ))
summary(thicDat2)



# Phon_Acc: Volume Stats ------------------------------------------------------------

# testing it with scaled Phon_Acc


# Check with fixed effects structure best explains variability in the data

m1 <- lmer(value ~ 
             1 
           + (1 | PID)
           , data = volDat)
m2 <- lmer(value ~ 
             Gyrus # Does including gyrus labels explain volumetric variability the than not including it?
           + (1 | PID)
           , data = volDat) 
m3 <- lmer(value ~ 
             Gyrus 
           + Hemisphere # Does including hemisphere labels better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m4 <- lmer(value ~ 
             Gyrus 
           + Hemisphere
           + Phon_Acc # Does including phonological accuracy better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m5 <- lmer(value ~ 
             Gyrus 
           * Hemisphere # Does an interaction between hemisphere and gyrus labels better explain the volumetric variability in the data?
           + Phon_Acc 
           + (1 | PID)
           , data = volDat) 
m6 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Phon_Acc 
           * Gyrus # Does an interaction between phonological accuracy and gyrus labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m7 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Phon_Acc 
           * Gyrus 
           + Phon_Acc 
           * Hemisphere # Does an interaction between phonological accuracy and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m8 <- lmer(value ~ 
             Gyrus 
           * Hemisphere
           * Phon_Acc # Does an interaction between phonological accuracy, gyrus, and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
anova(m1,m2,m3,m4,m5,m6,m7,m8)

#     Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# m1    3 6828967 6829000 -3414480  6828961                             
# m2   27 6385479 6385774 -3192713  6385425 4.4354e+05 24     <2e-16 ***
# m3   28 6381804 6382110 -3190874  6381748 3.6773e+03  1     <2e-16 ***
# m4   29 6381806 6382123 -3190874  6381748 1.0000e-04  1     0.9921    
# m5   52 6362358 6362926 -3181127  6362254 1.9494e+04 23     <2e-16 ***
# m6   76 6362029 6362859 -3180939  6361877 3.7689e+02 24     <2e-16 ***
# m7   77 6361963 6362804 -3180905  6361809 6.8001e+01  1     <2e-16 ***
# m8  100 6361861 6362953 -3180830  6361661 1.4844e+02 23     <2e-16 ***  

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m8), file = "volDat_Phon_Acc-scaled_model-summary.txt")

# Scaled is more interpretable, so we'll stick with that approach for the rest of the measures


# Manner_Acc: Volume Stats ------------------------------------------------------------

# Find best-fit model for Fixed Effects structure using PID as a random intercept to account for participant variability

m1 <- lmer(value ~ 
             1 
           + (1 | PID)
           , data = volDat)
m2 <- lmer(value ~ 
             Gyrus # Does including gyrus labels explain volumetric variability the than not including it?
           + (1 | PID)
           , data = volDat) 
m3 <- lmer(value ~ 
             Gyrus 
           + Hemisphere # Does including hemisphere labels better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m4 <- lmer(value ~ 
             Gyrus 
           + Hemisphere
           + Manner_Acc # Does including manner accuracy better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m5 <- lmer(value ~ 
             Gyrus 
           * Hemisphere # Does an interaction between hemisphere and gyrus labels better explain the volumetric variability in the data?
           + Manner_Acc 
           + (1 | PID)
           , data = volDat) 
m6 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Manner_Acc 
           * Gyrus # Does an interaction between manner accuracy and gyrus labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m7 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Manner_Acc 
           * Gyrus 
           + Manner_Acc 
           * Hemisphere # Does an interaction between manner accuracy and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m8 <- lmer(value ~ 
             Gyrus 
           * Hemisphere
           * Manner_Acc # Does an interaction between manner accuracy, gyrus, and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
anova(m1,m2,m3,m4,m5,m6,m7,m8)

#     Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# m1    3 6828967 6829000 -3414480  6828961                             
# m2   27 6385479 6385774 -3192713  6385425 443535.395 24  < 2.2e-16 ***
# m3   28 6381804 6382110 -3190874  6381748   3677.301  1  < 2.2e-16 ***
# m4   29 6381806 6382123 -3190874  6381748      0.000  1      0.995    
# m5   52 6362358 6362926 -3181127  6362254  19494.203 23  < 2.2e-16 ***
# m6   76 6362161 6362991 -3181004  6362009    245.060 24  < 2.2e-16 ***
# m7   77 6362118 6362959 -3180982  6361964     44.529  1  2.507e-11 ***
# m8  100 6362070 6363162 -3180935  6361870     93.999 23  1.501e-10 ***

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m8), file = "volDat_Manner_Acc-scaled_model-summary.txt")

# Scaled is more interpretable, so we'll stick with that approach for the rest of the measures


# Place_Acc: Volume Stats ------------------------------------------------------------

# Find best-fit model for Fixed Effects structure using PID as a random intercept to account for participant variability

m1 <- lmer(value ~ 
             1 
           + (1 | PID)
           , data = volDat)
m2 <- lmer(value ~ 
             Gyrus # Does including gyrus labels explain volumetric variability the than not including it?
           + (1 | PID)
           , data = volDat) 
m3 <- lmer(value ~ 
             Gyrus 
           + Hemisphere # Does including hemisphere labels better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m4 <- lmer(value ~ 
             Gyrus 
           + Hemisphere
           + Place_Acc # Does including Place accuracy better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m5 <- lmer(value ~ 
             Gyrus 
           * Hemisphere # Does an interaction between hemisphere and gyrus labels better explain the volumetric variability in the data?
           + Place_Acc 
           + (1 | PID)
           , data = volDat) 
m6 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Place_Acc 
           * Gyrus # Does an interaction between Place accuracy and gyrus labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m7 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Place_Acc 
           * Gyrus 
           + Place_Acc 
           * Hemisphere # Does an interaction between Place accuracy and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) # didn't converge
m8 <- lmer(value ~ 
             Gyrus 
           * Hemisphere
           * Place_Acc # Does an interaction between Place accuracy, gyrus, and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
anova(m1,m2,m3,m4,m5,m6,m7,m8)

#     Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# m1    3 6828967 6829000 -3414480  6828961                             
# m2   27 6385479 6385774 -3192713  6385425 443535.395 24  < 2.2e-16 ***
# m3   28 6381804 6382110 -3190874  6381748   3677.301  1  < 2.2e-16 ***
# m4   29 6381806 6382123 -3190874  6381748      0.000  1     0.9947    
# m5   52 6362358 6362926 -3181127  6362254  19494.203 23  < 2.2e-16 ***
# m6   76 6362164 6362994 -3181006  6362012    241.690 24  < 2.2e-16 ***
# m7   77 6362125 6362966 -3180985  6361971     41.579  1  1.132e-10 ***
# m8  100 6362076 6363168 -3180938  6361876     94.981 23  1.022e-10 ***

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m5), file = "volDat_Place_Acc-scaled_model-summary.txt")


# Voicing_Acc: Volume Stats ------------------------------------------------------------

# Find best-fit model for Fixed Effects structure using PID as a random intercept to account for participant variability

m1 <- lmer(value ~ 
             1 
           + (1 | PID)
           , data = volDat)
m2 <- lmer(value ~ 
             Gyrus # Does including gyrus labels explain volumetric variability the than not including it?
           + (1 | PID)
           , data = volDat) 
m3 <- lmer(value ~ 
             Gyrus 
           + Hemisphere # Does including hemisphere labels better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m4 <- lmer(value ~ 
             Gyrus 
           + Hemisphere
           + Voicing_Acc # Does including Voicing accuracy better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m5 <- lmer(value ~ 
             Gyrus 
           * Hemisphere # Does an interaction between hemisphere and gyrus labels better explain the volumetric variability in the data?
           + Voicing_Acc 
           + (1 | PID)
           , data = volDat) 
m6 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Voicing_Acc 
           * Gyrus # Does an interaction between Voicing accuracy and gyrus labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m7 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Voicing_Acc 
           * Gyrus 
           + Voicing_Acc 
           * Hemisphere # Does an interaction between Voicing accuracy and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m8 <- lmer(value ~ 
             Gyrus 
           * Hemisphere
           * Voicing_Acc # Does an interaction between Voicing accuracy, gyrus, and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
anova(m1,m2,m3,m4,m5,m6,m7,m8)

#     Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# m1    3 6828967 6829000 -3414480  6828961                             
# m2   27 6385479 6385774 -3192713  6385425 4.4354e+05 24  < 2.2e-16 ***
# m3   28 6381804 6382110 -3190874  6381748 3.6773e+03  1  < 2.2e-16 ***
# m4   29 6381806 6382123 -3190874  6381748 0.0000e+00  1    0.99577    
# m5   52 6362358 6362926 -3181127  6362254 1.9494e+04 23  < 2.2e-16 ***
# m6   76 6362303 6363133 -3181075  6362151 1.0336e+02 24  7.991e-12 ***
# m7   77 6362298 6363139 -3181072  6362144 6.2497e+00  1    0.01242 *  
# m8  100 6362318 6363410 -3181059  6362118 2.5974e+01 23    0.30210 

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m7), file = "volDat_Voicing_Acc-scaled_model-summary.txt")

