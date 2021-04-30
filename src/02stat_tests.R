colnames(rvm_gsh)
outliers15 <- which(rvm_gsh$salamander_id=='CONS3')
outliers15
outliers18 <- c(which(rvm_gsh$salamander_id=='CONS3'),which(rvm_gsh$salamander_id=='DS7'),which(rvm_gsh$salamander_id=='CHLS4'))
outliers18

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


# not an anova since no interaction, do ttests
which_controls <- which(rvm_gsh$treatment=='C')
which_24d <- which(rvm_gsh$treatment=='24D')
which_cps <- which(rvm_gsh$treatment=='CPS')
which15_controls <- which(rvm_gsh15$treatment=='C')
which15_24d <- which(rvm_gsh15$treatment=='24D')
which15_cps <- which(rvm_gsh15$treatment=='CPS')
which18_controls <- which(rvm_gsh18$treatment=='C')
which18_24d <- which(rvm_gsh18$treatment=='24D')
which18_cps <- which(rvm_gsh18$treatment=='CPS')
which_logged_controls <- which(rvm_gsh_logged$treatment=='C')
which_logged_24d <- which(rvm_gsh_logged$treatment=='24D')
which_logged_cps <- which(rvm_gsh_logged$treatment=='CPS')

###Welch Two Sample t-test
###1:5 24d
#1:5 with outliers
t.test(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d])
#data:  rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d]
#t = -1.0087, df = 19.697, p-value = 0.3254
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -7.542958  2.629019
#sample estimates:
#  mean of x mean of y 
#12.23892  14.69589 
#1:5 drop outliers
t.test(rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls],rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_24d])
#data:  rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls] and rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_24d]
#t = -2.0327, df = 17.403, p-value = 0.05763
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.1055902  0.1437023
#sample estimates:
#  mean of x mean of y 
#10.71495  14.69589 

###Welch Two Sample t-test
###1:8 24d
#1:8 with outliers
t.test(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d])
#data:  rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d]
#t = -0.92174, df = 17.605, p-value = 0.3691
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.691956  3.396786
#sample estimates:
#  mean of x mean of y 
#15.33464  17.98223 
#1:8 drop outliers
t.test(rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls],rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_24d])
#data:  rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls] and rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_24d]
#t = -3.4439, df = 17.898, p-value = 0.002915
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -9.573184 -2.316782
#sample estimates:
#  mean of x mean of y 
#13.25235  19.19733 

###Welch Two Sample t-test logged
#1:5 with outliers logged
t.test(log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d]))
#data:  log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d])
#t = -1.0713, df = 19.987, p-value = 0.2968
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.6010847  0.1931923
#sample estimates:
#  mean of x mean of y 
#2.408159  2.612106 
#1:8 with outliers logged
t.test(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d])
#data:  rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d]
#t = -0.92174, df = 17.605, p-value = 0.3691
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.691956  3.396786
#sample estimates:
#  mean of x mean of y 
#15.33464  17.98223 


#######Chlorpyrifos
###Welch Two Sample t-test
###1:5 chlorpyrifos
#with outliers
t.test(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps])
#data:  rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps]
#t = 0.81271, df = 17.337, p-value = 0.4274
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.832306  6.390034
#sample estimates:
#  mean of x mean of y 
#12.23892  10.46006

#drop outliers
t.test(rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls],rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_cps])
#data:  rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls] and rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_cps]
#t = 0.15531, df = 18.984, p-value = 0.8782
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3.180201  3.689980
#sample estimates:
#  mean of x mean of y 
#10.71495  10.46006 

###1:8
#with outliers
t.test(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps])
#data:  rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps]
#t = 0.62645, df = 16.122, p-value = 0.5398
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -4.107496  7.556417
#sample estimates:
#  mean of x mean of y 
#15.33464  13.61018 

#drop outliers
t.test(rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls],rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_cps])
#data:  rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls] and rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_cps]
#t = 0.42487, df = 17.162, p-value = 0.6762
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.727175  4.103762
#sample estimates:
#  mean of x mean of y 
#13.25235  12.56405


### LOGGED results
###Welch Two Sample t-test logged
###1:5 chlorpyrifos
#with or without outliers logged (outlier is 24d)
t.test(log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps]))
#data:  log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps])
#t = 0.69443, df = 18.839, p-value = 0.4959
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2415983  0.4813059
#sample estimates:
#  mean of x mean of y 
#2.408159  2.288306 
###1:8
#with or without outliers logged (outlier is 24d)
t.test(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps])
#data:  rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps]
#t = 0.62645, df = 16.122, p-value = 0.5398
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -4.107496  7.556417
#sample estimates:
#  mean of x mean of y 
#15.33464  13.61018 

#24d with outliers
t.test(log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d]))
#data:  log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d])
#t = -1.0713, df = 19.987, p-value = 0.2968
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.6010847  0.1931923
#sample estimates:
#  mean of x mean of y 
#2.408159  2.612106 
t.test(log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d]))
#data:  log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d])
#t = -1.0812, df = 19.895, p-value = 0.2925
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5678654  0.1802275
#sample estimates:
#  mean of x mean of y 
#2.635863  2.829682 

###
#24d with logged outlier dropped
t.test(log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_logged_controls]),log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_logged_24d]))
#data:  log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_logged_controls]) and log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_logged_24d])
#t = -1.163, df = 18.786, p-value = 0.2594
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.6458058  0.1846947
#sample estimates:
#  mean of x mean of y 
#2.408159  2.638715 

t.test(log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_logged_controls]),log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_logged_24d]))
#data:  log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_logged_controls]) and log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_logged_24d])
#t = -1.0699, df = 18.884, p-value = 0.2981
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5950722  0.1925897
#sample estimates:
#  mean of x mean of y 
#2.635863  2.837105 