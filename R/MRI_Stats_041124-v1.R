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

setwd('/Users/tifanibiro/Library/CloudStorage/GoogleDrive-tifbiro@gmail.com/My Drive/Research/Aphasia/lvPhon/Scripting/Brendan/Rebased/lvPhon/R/Resources')
data <- as.data.frame(read.csv(file='lvPhon-MRI-Acc-dat-v3.csv', header=TRUE, sep=','))
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

# First running it on Phon_Acc Grouped

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
           + Phon_Acc_Group # Does including phonological accuracy better explain volumetric variability than not including it?
           + (1 | PID)
           , data = volDat) 
m5 <- lmer(value ~ 
             Gyrus 
           * Hemisphere # Does an interaction between hemisphere and gyrus labels better explain the volumetric variability in the data?
           + Phon_Acc_Group 
           + (1 | PID)
           , data = volDat) 
m6 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Phon_Acc_Group 
           * Gyrus # Does an interaction between phonological accuracy and gyrus labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m7 <- lmer(value ~ 
             Gyrus 
           * Hemisphere 
           + Phon_Acc_Group 
           * Gyrus 
           + Phon_Acc_Group 
           * Hemisphere # Does an interaction between phonological accuracy and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
m8 <- lmer(value ~ 
             Gyrus 
           * Hemisphere
           * Phon_Acc_Group # Does an interaction between phonological accuracy, gyrus, and hemisphere labels better explain the volumetric variability in the data?
           + (1 | PID)
           , data = volDat) 
anova(m1,m2,m3,m4,m5,m6,m7,m8)

#     Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# m1   3 47483 47501 -23738    47477                               
# m2  27 44536 44696 -22241    44482 2995.542     24  < 2.2e-16 ***
# m3  28 44515 44682 -22230    44459   22.534      1  2.064e-06 ***
# m4  35 44510 44718 -22220    44440   19.436      7  0.0069245 ** 
# m5  58 44422 44767 -22153    44306  133.638     23  < 2.2e-16 ***
# m6 226 44573 45917 -22060    44121  185.058    168  0.1744061    
# m7 233 44560 45946 -22047    44094   26.991      7  0.0003346 *** # best fit model
# m8 394 44822 47166 -22017    44034   59.679    161  1.0000000   

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and phonological accuracy

# Write summary to a text file
capture.output(summary(m7), file = "volDat_Phon_Acc-ranked_model-summary.txt")

# Now testing it with scaled Phon_Acc


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
# m1   3 47483 47501 -23738    47477                                
# m2  27 44536 44696 -22241    44482 2995.5416     24  < 2.2e-16 *** # Main effect of gyrus
# m3  28 44515 44682 -22230    44459   22.5342      1  2.064e-06 *** # Main effect of hemisphere
# m4  29 44513 44685 -22227    44455    4.4608      1    0.03468 * # Main effect of phonological accuracy
# m5  52 44425 44734 -22160    44321  133.6380     23  < 2.2e-16 *** # Interaction between hemisphere and gyrus
# m6  76 44437 44889 -22142    44285   36.0760     24    0.05396 . # Could make an arguement for there being an interaction between phonological accuracy and gyrus, since it's p = .05
# m7  77 44436 44894 -22141    44282    2.5551      1    0.10994    
# m8 100 44473 45068 -22137    44273    9.1157     23    0.99555    

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m6), file = "volDat_Phon_Acc-scaled_model-summary.txt")

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
# m1   3 47483 47501 -23738    47477                                
# m2  27 44536 44696 -22241    44482 2995.5416     24  < 2.2e-16 ***
# m3  28 44515 44682 -22230    44459   22.5342      1  2.064e-06 ***
# m4  29 44515 44687 -22228    44457    2.4754      1     0.1156    
# m5  52 44427 44736 -22161    44323  133.6377     23  < 2.2e-16 ***
# m6  76 44447 44899 -22148    44295   27.9139     24     0.2637    
# m7  77 44446 44904 -22146    44292    2.5180      1     0.1126    
# m8 100 44485 45079 -22142    44285    7.9261     23     0.9985  

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m5), file = "volDat_Manner_Acc-scaled_model-summary.txt")

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
# m1   3 47483 47501 -23738    47477                                
# m2  27 44536 44696 -22241    44482 2995.5416     24  < 2.2e-16 ***
# m3  28 44515 44682 -22230    44459   22.5342      1  2.064e-06 ***
# m4  29 44514 44687 -22228    44456    2.6110      1     0.1061    
# m5  52 44427 44736 -22161    44323  133.6376     23  < 2.2e-16 ***
# m6  76 44445 44897 -22147    44293   29.5613     24     0.1997    
# m7  77 44445 44903 -22146    44291    1.9787      1     0.1595    
# m8 100 44483 45078 -22142    44283    7.8425     23     0.9986  

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
# m1   3 47483 47501 -23738    47477                                
# m2  27 44536 44696 -22241    44482 2995.5416     24  < 2.2e-16 ***
# m3  28 44515 44682 -22230    44459   22.5342      1  2.064e-06 ***
# m4  29 44511 44683 -22226    44453    6.2069      1   0.012725 *  
# m5  52 44423 44732 -22160    44319  133.6365     23  < 2.2e-16 ***
# m6  76 44421 44873 -22134    44269   50.3954     24   0.001261 ** 
# m7  77 44422 44880 -22134    44268    0.6756      1   0.411110    
# m8 100 44459 45054 -22130    44259    8.9367     23   0.996167 

# The best-fit model for explaining variability in the data includes an interaction between hemisphere and gyrus 

# Write summary to a text file
capture.output(summary(m6), file = "volDat_Voicing_Acc-scaled_model-summary.txt")

