# not an anova since no interaction, do ttests
# all the data
which_controls <- which(rvm_gsh$treatment=='C')
which_24d <- which(rvm_gsh$treatment=='24D')
which_cps <- which(rvm_gsh$treatment=='CPS')

which15_controls <- which(rvm_gsh15$treatment=='C')
which15_24d <- which(rvm_gsh15$treatment=='24D')
which15_cps <- which(rvm_gsh15$treatment=='CPS')

which18_controls <- which(rvm_gsh18$treatment=='C')
which18_24d <- which(rvm_gsh18$treatment=='24D')
which18_cps <- which(rvm_gsh18$treatment=='CPS')

# this is for logged data with outliers dropped (not just logged data)
# same for 1:5 and 1:8
which_logged_controls <- which(rvm_gsh_logged$treatment=='C')
which_logged_24d <- which(rvm_gsh_logged$treatment=='24D')
which_logged_cps <- which(rvm_gsh_logged$treatment=='CPS')


### We ultimately used the means of gsh 1:5 and 1:8
# but this section examines the 1:5 and 1:8 results separately

# exposures for 24d and chlorpyrifos
#24D significant for gsh when an outlier is dropped, for the 1:8 dilution but not 1:5
# chlorpyrifos never significant

#There are 16 total scenarios on this script:
# 2 dilutions * with or without outliers * logged/notlogged * 24d/chlorpyrifos
# we use none of them because we average the dilutions on the next script
# and use those for the manuscript


###Welch Two Sample t-test
### SCENARIO 1
###1:5 24d with outliers #NOT SIGNIFICANT
rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]
rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d]
t.test(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d])
#data:  rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d]
#t = -1.0087, df = 19.697, p-value = 0.3254
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -7.542958  2.629019
#sample estimates:
#  mean of x mean of y 
#12.23892  14.69589 

### SCENARIO 2
#1:5 24D drop outliers #ALMOST SIGNIFICANT
rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls]
rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_24d]
t.test(rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls],rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_24d])
#data:  rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls] and rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_24d]
#t = -2.0327, df = 17.403, p-value = 0.05763
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.1055902  0.1437023
#sample estimates:
#  mean of x mean of y 
#10.71495  14.69589 

### SCENARIO 3
###1:8 24d with outliers #NOT SIGNIFICANT
rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]
rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d]
t.test(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d])
#data:  rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d]
#t = -0.92174, df = 17.605, p-value = 0.3691
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.691956  3.396786
#sample estimates:
#  mean of x mean of y 
#15.33464  17.98223 

### SCENARIO 4
### 1:8 24D drop outliers #SIGNIFICANT
rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls]
rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_24d]
t.test(rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls],rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_24d])
#data:  rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls] and rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_24d]
#t = -3.4439, df = 17.898, p-value = 0.002915
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -9.573184 -2.316782
#sample estimates:
#  mean of x mean of y 
#13.25235  19.19733 

### SCENARIO 5
### 1:5 24D Logged with outliers # NOT SIGNIFICANT
log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls])
log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d])
t.test(log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d]))
#data:  log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_24d])
#t = -1.0713, df = 19.987, p-value = 0.2968
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.6010847  0.1931923
#sample estimates:
#  mean of x mean of y 
#2.408159  2.612106 

### SCENARIO 6
### 1:8 24D Logged with outliers logged # NOT SIGNIFICANT
log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls])
log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d])
t.test(log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d]))
#data:  log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_24d])
#t = -1.0812, df = 19.895, p-value = 0.2925
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5678654  0.1802275
#sample estimates:
#  mean of x mean of y 
#2.635863  2.829682  

### SCENARIO 7
# 1:5 24d with logged outlier dropped #ALMOST SIGNIFICANT
log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_controls])
log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_24d])
t.test(log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_controls]),log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_24d]))
#data:  log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_controls]) and log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_24d])
#t = -1.9157, df = 16.789, p-value = 0.07259
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.65288297  0.03179102
#sample estimates:
#  mean of x mean of y 
#2.408159  2.718705 

### SCENARIO 8
# 1:8 24d with logged outlier dropped # ALMOST SIGNIFICANT
log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_controls])
log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_24d])
t.test(log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_controls]),log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_24d]))
#data:  log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_controls]) and log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_24d])
#t = -2.0471, df = 14.611, p-value = 0.05907
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.61403264  0.01310285
#sample estimates:
#  mean of x mean of y 
#2.635863  2.936328  

#######Chlorpyrifos
### SCENARIO 9
###1:5 chlorpyrifos
# with outliers # NOT SIGNIFICANT
rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]
rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps]
t.test(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps])
#data:  rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps]
#t = 0.81271, df = 17.337, p-value = 0.4274
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.832306  6.390034
#sample estimates:
#  mean of x mean of y 
#12.23892  10.46006

### SCENARIO 10
# 1:5 chlorpyrifos drop outliers #NOT SIGNIFICANT
rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls]
rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_cps]
t.test(rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls],rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_cps])
#data:  rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_controls] and rvm_gsh15$gsh_1_5_dilution_nM_mL[which15_cps]
#t = 0.15531, df = 18.984, p-value = 0.8782
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3.180201  3.689980
#sample estimates:
#  mean of x mean of y 
#10.71495  10.46006 

### SCENARIO 11
### chlorpyrifos 1:8 with outliers #NOT SIGNIFICANT
rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]
rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps]
t.test(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls],rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps])
#data:  rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls] and rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps]
#t = 0.62645, df = 16.122, p-value = 0.5398
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -4.107496  7.556417
#sample estimates:
#  mean of x mean of y 
#15.33464  13.61018 

### SCENARIO 12
# chlorpyrifos 1:8 drop outliers # NOT SIGNIFICANT
rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls]
rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_cps]
t.test(rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls],rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_cps])
#data:  rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_controls] and rvm_gsh18$gsh_1_8_dilution_nM_mL[which18_cps]
#t = 0.42487, df = 17.162, p-value = 0.6762
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.727175  4.103762
#sample estimates:
#  mean of x mean of y 
#13.25235  12.56405

### SCENARIO 13
###1:5 chlorpyrifos logged with outliers # NOT SIGNIFICANT
log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls])
log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps])
t.test(log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps]))
#data:  log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_5_dilution_nM_mL[which_cps])
#t = 0.69443, df = 18.839, p-value = 0.4959
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2415983  0.4813059
#sample estimates:
#  mean of x mean of y 
#2.408159  2.288306 

### SCENARIO 14
###1:8 chlorpyrifos logged with outliers # NOT SIGNIFICANT
log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls])
log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps])
t.test(log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]),log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps]))
#data:  log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps])
#t = 0.45006, df = 18.499, p-value = 0.6579
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2697248  0.4171536
#sample estimates:
#  mean of x mean of y 
#2.635863  2.562149 

### SCENARIO 15
###1:5 chlorpyrifos logged drop outliers # NOT SIGNIFICANT (redundant, no outliers)
log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_controls])
log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_cps])
t.test(log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_controls]),log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_cps]))
#data:  log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_controls]) and log(rvm_gsh_logged$gsh_1_5_dilution_nM_mL[which_logged_cps])
#t = 0.69443, df = 18.839, p-value = 0.4959
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2415983  0.4813059
#sample estimates:
#  mean of x mean of y 
#2.408159  2.288306 

### SCENARIO 16
###1:8 chlorpyrifos logged drop outliers # NOT SIGNIFICANT (redundant, no outliers)
log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_controls])
log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_cps])
t.test(log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_controls]),log(rvm_gsh_logged$gsh_1_8_dilution_nM_mL[which_logged_cps]))
#data:  log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_controls]) and log(rvm_gsh$gsh_1_8_dilution_nM_mL[which_cps])
#t = 0.45006, df = 18.499, p-value = 0.6579
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2697248  0.4171536
#sample estimates:
#  mean of x mean of y 
#2.635863  2.562149 



