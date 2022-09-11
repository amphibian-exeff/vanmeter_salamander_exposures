############
# not an anova since no interaction, do ttests
# indices for the 3 scenarios: all data, dropping outlier, logged with dropped outlier
# all data
which_controls <- which(rvm_gsh$treatment=='C')
which_24d <- which(rvm_gsh$treatment=='24D')
which_cps <- which(rvm_gsh$treatment=='CPS')
# dropping outliers for averaged dilutions
which158_controls <- which(rvm_gsh_158$treatment=='C')
which158_24d <- which(rvm_gsh_158$treatment=='24D')
which158_cps <- which(rvm_gsh_158$treatment=='CPS')
# dropping logged outliers for averaged dilutions
which158_logged_controls <- which(rvm_gsh158_logged$treatment=='C')
which158_logged_24d <- which(rvm_gsh158_logged$treatment=='24D')
which158_logged_cps <- which(rvm_gsh158_logged$treatment=='CPS')


### We used the means of gsh 1:5 and 1:8 to generate the results
# this is what we use in the manuscript

# original data
View(rvm_gsh)
temp_view <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "C", "24D", "CPS"))

### SCENARIO 1
###averaged gsh 24d for 24d with outliers #NOT SIGNIFICANT
rvm_gsh$gsh_nM_mL[which_controls]
rvm_gsh$gsh_nM_mL[which_24d]
t.test(rvm_gsh$gsh_nM_mL[which_controls],rvm_gsh$gsh_nM_mL[which_24d])
#data:  rvm_gsh$gsh_nM_mL[which_controls] and rvm_gsh$gsh_nM_mL[which_24d]
#t = -0.9799, df = 18.409, p-value = 0.3398
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.015691  2.911137
#sample estimates:
#  mean of x mean of y 
#13.78678  16.33906 

### SCENARIO 1 #SIGNIFICANT
#averaged gsh for 24d  drop outliers
rvm_gsh_158$gsh_nM_mL[which158_controls] # high value of 31 dropped
rvm_gsh_158$gsh_nM_mL[which158_24d] # nothing dropped
t.test(rvm_gsh_158$gsh_nM_mL[which158_controls],rvm_gsh_158$gsh_nM_mL[which158_24d])
#data:  rvm_gsh_158$gsh_nM_mL[which158_controls] and rvm_gsh_158$gsh_nM_mL[which158_24d]
#t = -2.2871, df = 16.175, p-value = 0.03599
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -9.0154610 -0.3459656
#sample estimates:
#  mean of x mean of y 
#11.98365  16.66436 

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

###################
###################
#####################
# these are used for the manuscript

###gsh mean of 1_5 and 1_8

###Welch Two Sample t-test
###24d
#24d with outliers
rvm_gsh$gsh_nM_mL[which158_controls]
rvm_gsh$gsh_nM_mL[which158_24d]
t.test(rvm_gsh$gsh_nM_mL[which158_controls],rvm_gsh$gsh_nM_mL[which158_24d])
#data:  rvm_gsh$gsh_nM_mL[which158_controls] and rvm_gsh$gsh_nM_mL[which158_24d]
#t = -0.64131, df = 17.274, p-value = 0.5297
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -7.987325  4.260055
#sample estimates:
#  mean of x mean of y 
#14.26265  16.12628 

#24d drop outlier
rvm_gsh_158$gsh_nM_mL[which158_controls]
rvm_gsh_158$gsh_nM_mL[which158_24d]
t.test(rvm_gsh_158$gsh_nM_mL[which158_controls],rvm_gsh_158$gsh_nM_mL[which158_24d])
#data:  rvm_gsh_158$gsh_nM_mL[which158_controls] and rvm_gsh_158$gsh_nM_mL[which158_24d]
#t = -2.2871, df = 16.175, p-value = 0.03599
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -9.0154610 -0.3459656
#sample estimates:
#  mean of x mean of y 
#11.98365  16.66436  

#24d logged
t.test(log(rvm_gsh$gsh_nM_mL[which158_controls]),log(rvm_gsh$gsh_nM_mL[which_24d]))
#data:  log(rvm_gsh$gsh_nM_mL[which_controls]) and log(rvm_gsh$gsh_nM_mL[which_cps])
#t = 0.55897, df = 18.472, p-value = 0.5829
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2546375  0.4397163
#sample estimates:
#  mean of x mean of y 
#2.529329  2.436790 

#24d logged and drop outlier
t.test(log(rvm_gsh_158$gsh_nM_mL[which158_controls]),log(rvm_gsh_158$gsh_nM_mL[which158_24d]))
