#command + option + T = run entire section
#command + shift + R = add new section
#command + option + O = fold all sections
#command + option + shift + O = unfold all sections

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

# Define helper functions
## For mean centering
myCenter = function(x) x - mean(x, na.rm=T)

## For logging values
Logit = function(x) log(x/(1-x))


# Set up datasets & check distributions ---------------------------------------------------------

setwd('/Users/tifanibiro/Library/CloudStorage/GoogleDrive-tifbiro@gmail.com/My Drive/Research/Aphasia/lvPhon/Scripting/Tif/Resources')
expdat <- as.data.frame(read.csv(file='all_data_pvm_acc5.csv', header=TRUE, sep=','))
head(expdat)

#Factors 
expdat$Session_ID <- factor(expdat$Session_ID, 
                              levels = c("0","1","2"))
expdat$PID <- factor(expdat$PID)
expdat$wabaq_start <- factor(expdat$wabaq_start, 
                             levels = c('91-100', '81-90', '71-80', '61-70', '41-50'))
expdat$Target_Phoneme_ID <- factor(expdat$Target_Phoneme_ID)
expdat$Target <- factor(expdat$Target)
expdat$Word_Sess_Code <- factor(expdat$Word_Sess_Code)
expdat$Target_Word_Pos <- factor(expdat$Target_Word_Pos)
expdat$Target_Syll_Env <- factor(expdat$Target_Syll_Env)
expdat$Target_Con_Cluster <- factor(expdat$Target_Con_Cluster)
expdat$Syllable_NumID <- factor(expdat$Syllable_NumID)

#Check
summary(expdat)


# Statistics Set-up --------------------------------------------------------------

#Set-up

#Contrast building
levels(expdat$Session_ID)

#Make contrasts
Baseline_vs_12week <- c(2, -1, -1)
Active_vs_Sham <- c(0, 1, -1)

cHelmert2 <- cbind(Baseline_vs_12week, Active_vs_Sham)


# WAB Scores --------------------------------------------------------------

expdatW <- expdat

m1 <- lmer(
  wab1_aq ~ 
    Session_ID +
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdatW)

m2 <- lmer(
  wab1_aq ~ 
    Session_ID +
    (1 | PID) +
    (0 + Session_ID | PID), 
  contrasts = list(Session_ID = cHelmert2), data=expdatW) # Overfit

relLik(m1, m2)
# m1 wins by a lot

#    AIC(x)    AIC(y)      diff    relLik 
# 15509.07 -70776.20  86285.27       Inf 

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
expdatW <- romr.fnc(m1, expdatW, trim = 2.5)
expdatW$n.removed #0 outliers removed
expdatW$percent.removed #0% of the data
expdatW<-expdatW$data

summary(m1)

#                               Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                  6.958e+01  5.730e+00 7.000e+00   12.14 5.87e-06 ***
# Session_IDBaseline_vs_12week 4.685e-01  3.081e-02 3.283e+03   15.21  < 2e-16 ***
# Session_IDActive_vs_Sham     1.150e+00  5.539e-02 3.283e+03   20.77  < 2e-16 ***

m1a <- lmer(
  wab1_aq ~ 
    1 + 
    (1 | PID), 
  data=expdatW)

m1b <- lmer(
  wab1_aq ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdatW)

anova(m1a, m1b)

#    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# m1a  3 16087 16106 -8040.7    16081                             
# m1b  5 15505 15536 -7747.6    15495 586.18      2  < 2.2e-16 ***

summary(m1b)
#                              Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                  6.958e+01  5.730e+00 7.000e+00   12.14 5.87e-06 ***
# Session_IDBaseline_vs_12week 4.685e-01  3.081e-02 3.283e+03   15.21  < 2e-16 *** # Baseline better than post-stim
# Session_IDActive_vs_Sham     1.150e+00  5.539e-02 3.283e+03   20.77  < 2e-16 *** # Active better than sham stim
 
# The sample size is so small that we can't trust this effect

# WAB NWF Subscore  -------------------------------------------------------

# Naming and word finding

m1 <- lmer(
  wab1_nwf_total ~ 
    Session_ID +
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdatW)

m2 <- lmer(
  wab1_nwf_total ~ 
    Session_ID +
    (1 | PID) +
    (0 + Session_ID | PID), 
  contrasts = list(Session_ID = cHelmert2), data=expdatW) # Overfit

relLik(m1, m2)
# m2 wins, but is overfit, so we'll go with m1

#   AIC(x)     AIC(y)       diff     relLik 
# -1446.022 -82446.353  81000.330        Inf 

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
expdatW <- romr.fnc(m1, expdatW, trim = 2.5)
expdatW$n.removed #119 outliers removed
expdatW$percent.removed #3.6% of the data
expdatW<-expdatW$data

summary(m2)

#                               Estimate Std. Error        df t value Pr(>|t|)
# (Intercept)                  70.852796   0.290553  0.150450 243.855    0.345
# Session_IDBaseline_vs_12week  0.451585   0.061317  0.009565   7.365    0.953
# Session_IDActive_vs_Sham      0.220728   0.156084  1.192626   1.414    0.363

m2a <- lmer(
  wab1_nwf_total ~ 
    1 + 
    (1 | PID)  +
    (0 + Session_ID | PID), 
  data=expdatW)

m2b <- lmer(
  wab1_nwf_total ~ 
    Session_ID + 
    (1 | PID) +
    (0 + Session_ID | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdatW)

anova(m2a, m2b)



# Accuracy Set-up ---------------------------------------------------------

expdat2 <- as.data.frame(read.csv(file='AllAccScores_DamLev_080123.csv', header=TRUE, sep=','))
head(expdat2)

#Factors 
expdat2$Session_ID <- factor(expdat2$Session_ID, 
                            levels = c("0","1","2"))
expdat2$PID <- factor(expdat2$PID)
expdat2$wabaq_start <- factor(expdat2$wabaq_start, 
                             levels = c('91-100', '81-90', '71-80', '61-70', '41-50'))
expdat2$Target_Phoneme_ID <- factor(expdat2$Target_Phoneme_ID)
expdat2$Target <- factor(expdat2$Target)
expdat2$Word_Sess_Code <- factor(expdat2$Word_Sess_Code)
expdat2$Target_Word_Pos <- factor(expdat2$Target_Word_Pos)
expdat2$Target_Syll_Env <- factor(expdat2$Target_Syll_Env)
expdat2$Target_Con_Cluster <- factor(expdat2$Target_Con_Cluster)
expdat2$Syllable_NumID <- factor(expdat2$Syllable_NumID)

#Check
summary(expdat2)


####### BROKEN ######
# Phon_Acc ----------------------------------------------------------------

expdat2PA <- expdat2

m1 <- lmer(
  Phon_Acc ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m2 <- lmer(
  Phon_Acc ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID) +
    (1 | Target), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2PA) # Overfit

relLik(m1, m2)
# m2 wins, but is overfit, so we'll go with m1 (also for consistency's sake, we should use this model)
#       AIC(x)      AIC(y)        diff      relLik 
# 4.039195e+03 3.901790e+03 1.374046e+02 6.871011e+29

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
expdat2PA <- romr.fnc(m1, expdat2PA, trim = 2.5)
expdat2PA$n.removed #0 outliers removed
expdat2PA$percent.removed #0% of the data
expdat2PA<-expdat2PA$data

summary(m1)
#                                                   Estimate Std. Error         df t value Pr(>|t|)   
# (Intercept)                                      9.507e-01  9.905e-02  3.165e+00   9.598  0.00192 **
# Session_IDBaseline_vs_12week                   1.200e-02  1.981e-02  3.272e+03   0.606  0.54452   
# Session_IDActive_vs_Sham                      -4.867e-02  3.304e-02  3.272e+03  -1.473  0.14081   
# wabaq_start81-90                                -5.205e-02  1.398e-01  3.138e+00  -0.372  0.73337   
# wabaq_start71-80                                -2.082e-01  1.211e-01  3.145e+00  -1.719  0.17988   
# wabaq_start61-70                                -2.725e-01  1.138e-01  3.105e+00  -2.394  0.09351 . 
# wabaq_start41-50                                -3.440e-01  1.391e-01  3.076e+00  -2.473  0.08770 . 
# Session_IDBaseline_vs_12week:wabaq_start81-90 -1.260e-02  2.768e-02  3.272e+03  -0.455  0.64898   
# Session_IDActive_vs_Sham:wabaq_start81-90      6.153e-02  4.456e-02  3.272e+03   1.381  0.16744   
# Session_IDBaseline_vs_12week:wabaq_start71-80 -2.725e-02  2.336e-02  3.272e+03  -1.167  0.24346   
# Session_IDActive_vs_Sham:wabaq_start71-80      9.339e-02  4.034e-02  3.272e+03   2.315  0.02065 * 
# Session_IDBaseline_vs_12week:wabaq_start61-70 -2.595e-02  2.112e-02  3.273e+03  -1.229  0.21922   
# Session_IDActive_vs_Sham:wabaq_start61-70      1.108e-01  3.583e-02  3.272e+03   3.093  0.00200 **
# Session_IDBaseline_vs_12week:wabaq_start41-50  8.193e-03  2.543e-02  3.272e+03   0.322  0.74730   
# Session_IDActive_vs_Sham:wabaq_start41-50      5.044e-02  4.195e-02  3.272e+03   1.202  0.22930  

m_Baseline <- lmer(
  Phon_Acc ~ 
    1 + 
    (1 | PID), 
  data=expdat2PA)

m_Sess <- lmer(
  Phon_Acc ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_WABS <- lmer(
  Phon_Acc ~ 
    Session_ID + 
    wabaq_start + 
    (1 | PID),
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_SessxWABS <- lmer(
  Phon_Acc ~ 
    Session_ID *
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2PA)

anova(m_Baseline, m_Sess, m_WABS, m_SessxWABS)
#     Df   AIC   BIC  logLik deviance     Chisq Chi Df Pr(>Chisq)    
# m_Baseline   3 3974.5 3992.8 -1984.3   3968.5                             
# m_Sess       5 3964.1 3994.6 -1977.0   3954.1 14.442      2   0.000731 ***
# m_WABS       9 3960.7 4015.5 -1971.3   3942.7 11.421      4   0.022214 *  
# m_SessxWABS 17 3959.0 4062.7 -1962.5   3925.0 17.607      8   0.024376 *   

summary(m_SessxWABS)

#same as above summary

# FeatureWeighted_PhonAcc ----------------------------------------------------------------

expdat2FWPA <- expdat2

m1 <- lmer(
  FeatureWeighted_PhonAcc ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2FWPA)

m2 <- lmer(
  FeatureWeighted_PhonAcc ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID) +
    (1 | Target), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2FWPA) # Overfit

relLik(m1, m2)

#       AIC(x)      AIC(y)        diff      relLik 
# -2.360637e+02 -4.487304e+02  2.126667e+02  1.513508e+46 

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
expdat2FWPA <- romr.fnc(m1, expdat2FWPA, trim = 2.5)
expdat2FWPA$n.removed #202 outliers removed
expdat2FWPA$percent.removed #6.14% of the data
expdat2FWPA<-expdat2FWPA$data

summary(m1)
#                                                   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                      9.725e-01  3.718e-02  3.291e+00  26.157 6.22e-05 ***
# Session_IDBaseline_vs_12week                   7.721e-03  1.031e-02  3.272e+03   0.748  0.45422    
# Session_IDActive_vs_Sham                      -2.506e-02  1.721e-02  3.272e+03  -1.457  0.14534    
# wabaq_start81-90                                -1.759e-02  5.236e-02  3.237e+00  -0.336  0.75758    
# wabaq_start71-80                                -7.771e-02  4.540e-02  3.251e+00  -1.712  0.17834    
# wabaq_start61-70                                -8.954e-02  4.254e-02  3.170e+00  -2.105  0.12108    
# wabaq_start41-50                                -1.387e-01  5.186e-02  3.113e+00  -2.675  0.07244 .  
# Session_IDBaseline_vs_12week:wabaq_start81-90 -1.196e-02  1.442e-02  3.272e+03  -0.829  0.40693    
# Session_IDActive_vs_Sham:wabaq_start81-90      4.390e-02  2.321e-02  3.272e+03   1.892  0.05863 .  
# Session_IDBaseline_vs_12week:wabaq_start71-80 -1.835e-02  1.217e-02  3.272e+03  -1.508  0.13162    
# Session_IDActive_vs_Sham:wabaq_start71-80      4.712e-02  2.101e-02  3.272e+03   2.243  0.02496 *  
# Session_IDBaseline_vs_12week:wabaq_start61-70 -1.798e-02  1.100e-02  3.274e+03  -1.635  0.10222    
# Session_IDActive_vs_Sham:wabaq_start61-70      4.755e-02  1.866e-02  3.272e+03   2.548  0.01087 *  
# Session_IDBaseline_vs_12week:wabaq_start41-50 -9.817e-03  1.324e-02  3.272e+03  -0.741  0.45853    
# Session_IDActive_vs_Sham:wabaq_start41-50      5.701e-02  2.185e-02  3.272e+03   2.610  0.00911 ** 

m_Baseline <- lmer(
  FeatureWeighted_PhonAcc ~ 
    1 + 
    (1 | PID), 
  data=expdat2PA)

m_Sess <- lmer(
  FeatureWeighted_PhonAcc ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_WABS <- lmer(
  FeatureWeighted_PhonAcc ~ 
    Session_ID + 
    wabaq_start + 
    (1 | PID),
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_SessxWABS <- lmer(
  FeatureWeighted_PhonAcc ~ 
    Session_ID *
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2PA)

anova(m_Baseline, m_Sess, m_WABS, m_SessxWABS)

#             Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# m_Baseline   3 -324.13 -305.84 165.07  -330.13                             
# m_Sess       5 -339.68 -309.18 174.84  -349.68 19.544      2  5.703e-05 ***
# m_WABS       9 -343.17 -288.28 180.59  -361.17 11.495      4    0.02153 *  
# m_SessxWABS 17 -339.14 -235.46 186.57  -373.14 11.964      8    0.15281    

summary(m_WABS)

#                                 Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                     9.677e-01  3.746e-02  3.283e+00  25.835  6.6e-05 ***
# Session_IDBaseline_vs_12week -7.237e-03  2.808e-03  3.282e+03  -2.577 0.010016 *  
# Session_IDActive_vs_Sham      1.898e-02  5.053e-03  3.281e+03   3.756 0.000176 ***
# wabaq_start81-90               -1.343e-02  5.271e-02  3.217e+00  -0.255 0.814245    
# wabaq_start71-80               -7.358e-02  4.574e-02  3.245e+00  -1.608 0.199249    
# wabaq_start61-70               -8.546e-02  4.287e-02  3.167e+00  -1.993 0.135363    
# wabaq_start41-50               -1.367e-01  5.223e-02  3.104e+00  -2.616 0.076507 .  

# PVMWeighted_PhonAcc ----------------------------------------------------------------

expdat2PhPA <- expdat2

m1 <- lmer(
  PVMWeighted_PhonAcc ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PhPA)

m2 <- lmer(
  PVMWeighted_PhonAcc ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID) +
    (1 | Target), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2PhPA)

relLik(m1, m2)

#        AIC(x)      AIC(y)        diff      relLik 
# 4.598197e+02 2.414632e+02 2.183565e+02 2.603279e+47 

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
expdat2PhPA <- romr.fnc(m1, expdat2PhPA, trim = 2.5)
expdat2PhPA$n.removed #170 outliers removed
expdat2PhPA$percent.removed #5.17% of the data
expdat2PhPA<-expdat2PhPA$data

summary(m1)
#                                                   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                      9.711e-01  5.139e-02  3.181e+00  18.898 0.000224 ***
# Session_IDBaseline_vs_12week                   9.017e-03  1.147e-02  3.272e+03   0.786 0.431789    
# Session_IDActive_vs_Sham                      -3.067e-02  1.913e-02  3.272e+03  -1.603 0.109012    
# wabaq_start81-90                                -2.584e-02  7.248e-02  3.147e+00  -0.357 0.744036    
# wabaq_start71-80                                -1.053e-01  6.281e-02  3.156e+00  -1.677 0.187669    
# wabaq_start61-70                                -1.192e-01  5.898e-02  3.106e+00  -2.021 0.133423    
# wabaq_start41-50                                -1.685e-01  7.202e-02  3.069e+00  -2.340 0.099266 .  
# Session_IDBaseline_vs_12week:wabaq_start81-90 -1.371e-02  1.603e-02  3.272e+03  -0.855 0.392508    
# Session_IDActive_vs_Sham:wabaq_start81-90      4.773e-02  2.581e-02  3.272e+03   1.850 0.064424 .  
# Session_IDBaseline_vs_12week:wabaq_start71-80 -1.913e-02  1.353e-02  3.272e+03  -1.414 0.157411    
# Session_IDActive_vs_Sham:wabaq_start71-80      5.093e-02  2.336e-02  3.272e+03   2.181 0.029276 *  
# Session_IDBaseline_vs_12week:wabaq_start61-70 -1.869e-02  1.223e-02  3.273e+03  -1.529 0.126481    
# Session_IDActive_vs_Sham:wabaq_start61-70      5.979e-02  2.074e-02  3.272e+03   2.882 0.003974 ** 
# Session_IDBaseline_vs_12week:wabaq_start41-50 -6.951e-03  1.472e-02  3.272e+03  -0.472 0.636868    
# Session_IDActive_vs_Sham:wabaq_start41-50      5.368e-02  2.429e-02  3.272e+03   2.210 0.027193 *  

m_Baseline <- lmer(
  PVMWeighted_PhonAcc ~ 
    1 + 
    (1 | PID), 
  data=expdat2PA)

m_Sess <- lmer(
  PVMWeighted_PhonAcc ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_WABS <- lmer(
  PVMWeighted_PhonAcc ~ 
    Session_ID + 
    wabaq_start + 
    (1 | PID),
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_SessxWABS <- lmer(
  PVMWeighted_PhonAcc ~ 
    Session_ID *
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2PA)

anova(m_Baseline, m_Sess, m_WABS, m_SessxWABS)

#             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# m_Baseline   3 372.53 390.83 -183.27   366.53                             
# m_Sess       5 360.76 391.25 -175.38   350.76 15.777      2  0.0003751 ***
# m_WABS       9 358.57 413.46 -170.28   340.57 10.187      4  0.0374001 *  
# m_SessxWABS 17 362.14 465.82 -164.07   328.14 12.430      8  0.1330371    

summary(m_WABS)
#                                  Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                     9.658e-01  5.192e-02  3.178e+00  18.600 0.000237 ***
# Session_IDBaseline_vs_12week -6.197e-03  3.124e-03  3.283e+03  -1.984 0.047354 *  
# Session_IDActive_vs_Sham      2.000e-02  5.618e-03  3.281e+03   3.560 0.000376 ***
# wabaq_start81-90               -2.041e-02  7.319e-02  3.137e+00  -0.279 0.797721    
# wabaq_start71-80               -1.006e-01  6.347e-02  3.154e+00  -1.585 0.206787    
# wabaq_start61-70               -1.151e-01  5.962e-02  3.106e+00  -1.931 0.145926    
# wabaq_start41-50               -1.649e-01  7.277e-02  3.066e+00  -2.266 0.106458  

# Damerau_Levenshtein ----------------------------------------------------------------

expdat2DL <- expdat2

m1 <- lmer(
  Damerau_Levenshtein ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2DL)

m2 <- lmer(
  Damerau_Levenshtein ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID) +
    (1 | Target), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2DL)

relLik(m1, m2)

#       AIC(x)      AIC(y)        diff      relLik 
#  1.377025e+04  1.291603e+04  8.542179e+02 3.097816e+185 

#check model assumptions
mcp.fnc(m2) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
expdat2DL <- romr.fnc(m2, expdat2DL, trim = 2.5)
expdat2DL$n.removed #90 outliers removed
expdat2DL$percent.removed #2.74% of the data
expdat2DL<-expdat2DL$data

summary(m2)
#                                                   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                      5.937e-02  8.478e-01  3.440e+00   0.070   0.9480  
# Session_IDBaseline_vs_12week                   4.554e-02  7.593e-02  3.253e+03   0.600   0.5487  
# Session_IDActive_vs_Sham                      -1.867e-02  1.268e-01  3.253e+03  -0.147   0.8829  
# wabaq_start81-90                                 1.250e-01  1.161e+00  3.028e+00   0.108   0.9210  
# wabaq_start71-80                                 8.828e-01  1.006e+00  3.030e+00   0.878   0.4441  
# wabaq_start61-70                                 2.115e+00  9.476e-01  3.022e+00   2.232   0.1111  
# wabaq_start41-50                                 2.173e+00  1.160e+00  3.017e+00   1.873   0.1573  
# Session_IDBaseline_vs_12week:wabaq_start81-90 -2.628e-01  1.063e-01  3.254e+03  -2.473   0.0134 *
# Session_IDActive_vs_Sham:wabaq_start81-90      4.196e-02  1.716e-01  3.254e+03   0.245   0.8068  
# Session_IDBaseline_vs_12week:wabaq_start71-80  4.833e-02  8.969e-02  3.254e+03   0.539   0.5900  
# Session_IDActive_vs_Sham:wabaq_start71-80     -1.422e-01  1.550e-01  3.254e+03  -0.917   0.3591  
# Session_IDBaseline_vs_12week:wabaq_start61-70  4.017e-02  8.122e-02  3.254e+03   0.495   0.6209  
# Session_IDActive_vs_Sham:wabaq_start61-70     -1.951e-01  1.375e-01  3.254e+03  -1.419   0.1561  
# Session_IDBaseline_vs_12week:wabaq_start41-50 -1.969e-01  9.772e-02  3.254e+03  -2.015   0.0440 *
# Session_IDActive_vs_Sham:wabaq_start41-50      1.014e-03  1.616e-01  3.255e+03   0.006   0.9950  

m_Baseline <- lmer(
  Damerau_Levenshtein ~ 
    1 + 
    (1 | PID) +
    (1 | Target), 
  data=expdat2PA)

m_Sess <- lmer(
  Damerau_Levenshtein ~ 
    Session_ID + 
    (1 | PID) +
    (1 | Target), 
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_WABS <- lmer(
  Damerau_Levenshtein ~ 
    Session_ID + 
    wabaq_start + 
    (1 | PID) +
    (1 | Target),
  contrasts = list(Session_ID = cHelmert2), 
  data=expdat2PA)

m_SessxWABS <- lmer(
  Damerau_Levenshtein ~ 
    Session_ID *
    wabaq_start + 
    (1 | PID) +
    (1 | Target), 
  contrasts = list(Session_ID = cHelmert2), data=expdat2PA)

anova(m_Baseline, m_Sess, m_WABS, m_SessxWABS)

#             Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# m_Baseline   4 12912 12937 -6452.1    12904                             
# m_Sess       6 12901 12938 -6444.6    12889 14.883      2  0.0005863 ***
# m_WABS      10 12899 12960 -6439.4    12879 10.424      4  0.0338665 *  
# m_SessxWABS 18 12885 12995 -6424.5    12849 29.958      8  0.0002150 ***

summary(m_SessxWABS)
# See above

summary(m_Sess)

# Silence and pause segments set up ----------------------------------------------

# look at pause rate (% of pauses = (pause duration/total duration)) improvement
# look at WAB spontaneous speech (SS) subscore improvement
# use another column as dependent variable?


pauses_data <- as.data.frame(read.csv(file='pauses.csv', header=TRUE, sep=','))

#pauses_data$PID <- as.integer(substr(pauses_data$PID, nchar(pauses_data$PID)-1, nchar(pauses_data$PID)))
# Above code brought up an error, so I commented it out. You actually can leave the characters in to use it as a random effect

# Factors

pauses_data$Session_ID <- factor(pauses_data$Session_ID, levels = c("Baseline","Active", "Sham"))
pauses_data$PID <- factor(pauses_data$PID)


#Contrast building
levels(pauses_data$Session_ID)

#Make contrasts
Baseline_vs_12week <- c(2, -1, -1)
Active_vs_Sham <- c(0, 1, -1)

cHelmert2 <- cbind(Baseline_vs_12week, Active_vs_Sham)


# WAB SS subscore ---------------------------------------------------------

m1 <- lmer(
  wab1_ss_total ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_data)

m2 <- lmer(
  wab1_ss_total ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID) +
    (0 + Session_ID | PID), 
  contrasts = list(Session_ID = cHelmert2), data=pauses_data) # Overfit

relLik(m1, m2)
# m2 wins, but is overfit, so we'll go with m1
#   AIC(x)     AIC(y)       diff     relLik 
# 8246.283 -90827.041  99073.324        Inf 

#check model assumptions
mcp.fnc(m1) 


#Remove outliers that are 2.5 standard deviations from 0 -- none
pauses_data <- romr.fnc(m1, pauses_data, trim = 2.5)
pauses_data$n.removed #122 outliers removed
pauses_data$percent.removed #3.25% of the data
pauses_data<-pauses_data$data

summary(m1)
#                                                                 Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)                                                      15.00260    2.68592    7.99879   5.586 0.000519**
# Session_IDBaseline_vs_12week                                   -0.27406    0.01655 3739.00124 -16.559  < 2e-16***
# Session_IDActive_vs_Sham                                       -0.28093    0.03034 3739.00164  -9.260  < 2e-16***
# wabaq_startImproved_Treat                                   3.33074    4.65222    7.99929   0.716 0.494391
# wabaq_startNo_Improved                                     -1.29538    3.00297    7.99903  -0.431 0.677592
# Session_IDBaseline_vs_12week:wabaq_startImproved_Treat    0.10740    0.03665 3739.00014   2.931 0.003405**
# Session_IDActive_vs_Sham:wabaq_startImproved_Treat        0.78093    0.05740 3739.00035  13.605  < 2e-16***
# Session_IDBaseline_vs_12week:wabaq_startNo_Improved       0.54921    0.02030 3739.02389  27.054  < 2e-16***
# Session_IDActive_vs_Sham:wabaq_startNo_Improved           0.58743    0.03441 3739.01124  17.071  < 2e-16***

m1a <- lmer(
  wab1_ss_total ~ 
    1 + 
    (1 | PID), 
  data=pauses_data)

m1b <- lmer(
  wab1_ss_total ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_data)

m1c <- lmer(
  wab1_ss_total ~ 
    Session_ID + 
    wabaq_start + 
    (1 | PID),
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_data)

m1d<- lmer(
  wab1_ss_total ~ 
    Session_ID *
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), data=pauses_data)

anova(m1a, m1b, m1c, m1d)
#     Df    AIC    BIC  logLik deviance     Chisq Chi Df Pr(>Chisq)    
# m1a  3 8762.1 8780.7 -4378.1   8756.1                                
# m1b  5 8341.9 8372.9 -4165.9   8331.9  424.2445      2     <2e-16 ***#Main effect of Session Type
# m1c  7 8343.8 8387.1 -4164.9   8329.8    2.1225      2      0.346    
# m1d 11 6134.6 6202.8 -3056.3   6112.6 2217.1630      4     <2e-16 ***#Interaction between Session Type and Improvement Group

# Pause Rate --------------------------------------------------------------

pauses_dataPR <- pauses_data
pauses_dataPR$wab1_aq <- factor(pauses_dataPR$wab1_aq)
pauses_dataPR$wab1_nwf_total <- factor(pauses_dataPR$wab1_nwf_total)
pauses_dataPR$wabaq_start <- factor(pauses_dataPR$wabaq_start)
pauses_dataPR$SS_wabaq_start <- factor(pauses_dataPR$SS_wabaq_start)


#Build models

m1 <- lmer(
  pause_rate ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_dataPR)

m2 <- lmer(
  pause_rate ~ 
    Session_ID * 
    wabaq_start + 
    (1 | PID) +
    (0 + Session_ID | PID), 
  contrasts = list(Session_ID = cHelmert2), data=pauses_dataPR) # Overfit

relLik(m1, m2)
# m2 wins with lowest AIC, but is overfit, so we'll go with m1
#    AIC(x)     AIC(y)       diff     relLik 
# -12037.39 -105248.62   93211.23        Inf 

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
pauses_dataPR <- romr.fnc(m1, pauses_dataPR, trim = 2.5)
pauses_dataPR$n.removed #100 outliers removed
pauses_dataPR$percent.removed #2.75% of the data
pauses_dataPR<-pauses_dataPR$data

summary(m1)
#                                                                  Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)                                                     6.196e-01  3.940e-02  7.984e+00  15.725 2.73e-07***
# Session_IDBaseline_vs_12week                                  2.218e-02  1.050e-03  3.617e+03  21.133  < 2e-16***
# Session_IDActive_vs_Sham                                      2.954e-02  1.924e-03  3.617e+03  15.356  < 2e-16***
# wabaq_startImproved_Treat                                -7.924e-02  6.827e-02  7.993e+00  -1.161  0.27921
# wabaq_startNo_Improved                                   -3.193e-02  4.406e-02  7.990e+00  -0.725  0.48937
# Session_IDBaseline_vs_12week:wabaq_startImproved_Treat -1.925e-03  2.324e-03  3.617e+03  -0.828  0.40745
# Session_IDActive_vs_Sham:wabaq_startImproved_Treat     -4.615e-04  3.640e-03  3.617e+03  -0.127  0.89913
# Session_IDBaseline_vs_12week:wabaq_startNo_Improved    -1.608e-02  1.317e-03  3.618e+03 -12.216  < 2e-16***
# Session_IDActive_vs_Sham:wabaq_startNo_Improved        -5.721e-03  2.200e-03  3.618e+03  -2.600  0.00936**

m1a <- lmer(
  pause_rate ~ 
    1 + 
    (1 | PID), 
  data=pauses_dataPR)

m1b <- lmer(
  pause_rate ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_dataPR)

m1c <- lmer(
  pause_rate ~ 
    Session_ID + 
    wabaq_start + 
    (1 | PID),
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_dataPR)

m1d<- lmer(
  pause_rate ~ 
    Session_ID *
    wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), data=pauses_dataPR)

anova(m1a, m1b, m1c, m1d)
#     Df    AIC    BIC logLik deviance     Chisq Chi Df Pr(>Chisq)    
# m1a  3 -11472 -11454 5739.3   -11478                                
# m1b  5 -12502 -12471 6255.9   -12512 1033.3148      2     <2e-16 *** #Main effect of session type
# m1c  7 -12499 -12456 6256.7   -12513    1.5954      2     0.4504    
# m1d 11 -12916 -12848 6468.8   -12938  424.1987      4     <2e-16 *** #Interation between session type and response group


# Pause Rate: Subscore group ----------------------------------------------------------

#Build models

m1 <- lmer(
  pause_rate ~ 
    Session_ID * 
    SS_wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_dataPR)

m2 <- lmer(
  pause_rate ~ 
    Session_ID * 
    SS_wabaq_start + 
    (1 | PID) +
    (0 + Session_ID | PID), 
  contrasts = list(Session_ID = cHelmert2), data=pauses_dataPR) # Overfit

relLik(m1, m2)
# m2 wins with lowest AIC, but is overfit, we'll go with m1 for consistency's sake
# AIC(x)     AIC(y)       diff     relLik 
# -12553.93 -106872.21   94318.28        Inf 

#check model assumptions
mcp.fnc(m1) 

#Remove outliers that are 2.5 standard deviations from 0 -- none
pauses_dataPR <- romr.fnc(m1, pauses_dataPR, trim = 2.5)
pauses_dataPR$n.removed #98 outliers removed
pauses_dataPR$percent.removed #2.77% of the data
pauses_dataPR<-pauses_dataPR$data

summary(m1)
#                                                                        Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                                           5.991e-01  3.782e-02  8.006e+00  15.843 2.50e-07 ***
# Session_IDBaseline_vs_12week                                        1.676e-02  1.234e-03  3.519e+03  13.585  < 2e-16 ***
# Session_IDActive_vs_Sham                                            3.199e-02  1.672e-03  3.517e+03  19.129  < 2e-16 ***
# SS_wabaq_startSS_Improved_Treat                                -7.251e-02  4.883e-02  8.013e+00  -1.485 0.175824    
# SS_wabaq_startSS_No_Improved                                    6.841e-03  4.367e-02  8.009e+00   0.157 0.879391    
# Session_IDBaseline_vs_12week:SS_wabaq_startSS_Improved_Treat -1.189e-02  1.892e-03  3.523e+03  -6.283 3.72e-10 ***
# Session_IDActive_vs_Sham:SS_wabaq_startSS_Improved_Treat     -2.432e-02  2.576e-03  3.518e+03  -9.441  < 2e-16 ***
# Session_IDBaseline_vs_12week:SS_wabaq_startSS_No_Improved    -7.070e-03  1.421e-03  3.519e+03  -4.976 6.82e-07 ***
# Session_IDActive_vs_Sham:SS_wabaq_startSS_No_Improved        -7.285e-03  1.962e-03  3.518e+03  -3.712 0.000208 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

m1a <- lmer(
  pause_rate ~ 
    1 + 
    (1 | PID), 
  data=pauses_dataPR)

m1b <- lmer(
  pause_rate ~ 
    Session_ID + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_dataPR)

m1c <- lmer(
  pause_rate ~ 
    Session_ID + 
    SS_wabaq_start + 
    (1 | PID),
  contrasts = list(Session_ID = cHelmert2), 
  data=pauses_dataPR)

m1d<- lmer(
  pause_rate ~ 
    Session_ID *
    SS_wabaq_start + 
    (1 | PID), 
  contrasts = list(Session_ID = cHelmert2), data=pauses_dataPR)

anova(m1a, m1b, m1c, m1d)
#     Df    AIC    BIC logLik deviance     Chisq Chi Df Pr(>Chisq)    
# m1a  3 -11386 -11367 5695.9   -11392                                
# m1b  5 -13205 -13174 6607.4   -13215 1822.9695      2    < 2e-16 *** #Main effect of session type
# m1c  7 -13207 -13164 6610.6   -13221    6.4425      2    0.03991 *  #Main effect of improvement group
# m1d 11 -13527 -13460 6774.5   -13549  327.8218      4    < 2e-16 *** #Interaction between session type and improvement group

#### NOTES ###
# random effects on right side of |
#0+Condition | Subject_ID)           
#(1 | Subject_ID)+(0+Condition | Subject_ID)
#(1 | Subject_ID)+(0+Condition | Subject_ID) = (1+Condition | Subject_ID)

