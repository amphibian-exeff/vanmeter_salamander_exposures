colnames(rvm_gsh)

#svl ~ body weight
svl_mm_model_wbw <- lm(svl_mm ~ weight_g, data=rvm_gsh)
summary(svl_mm_model_wbw)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  62.9125     4.0721  15.450 4.12e-16 ***
#  weight_g      1.2911     0.2177   5.931 1.49e-06 ***

#snout vent length
svl_mm_model <- lm(svl_mm ~ treatment, data=rvm_gsh)
summary(svl_mm_model)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    88.189      1.674  52.667   <2e-16 ***
#  treatmentC     -1.075      2.368  -0.454    0.653    
#treatmentCPS   -3.294      2.368  -1.391    0.174 

#body weight
bw_g_model <- lm(weight_g ~ treatment, data=rvm_gsh)
summary(bw_g_model)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   19.0275     0.9565  19.894   <2e-16 ***
#  treatmentC    -0.2839     1.3526  -0.210    0.835    
#treatmentCPS  -1.4486     1.3526  -1.071    0.293


