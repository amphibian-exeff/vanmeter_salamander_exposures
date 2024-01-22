# not an anova since no interaction, do ttests
which_controls <- which(rvm_gsh$treatment=='C')
which_24d <- which(rvm_gsh$treatment=='24D')
which_cps <- which(rvm_gsh$treatment=='CPS')
which_controls_ache <- which(rvm_gsh_ache$treatment=='C')
which_24d_ache <- which(rvm_gsh_ache$treatment=='24D')
which_cps_ache <- which(rvm_gsh_ache$treatment=='CPS')


which_logged_controls <- which(rvm_gsh_logged$treatment=='C')
which_logged_24d <- which(rvm_gsh_logged$treatment=='24D')
which_logged_cps <- which(rvm_gsh_logged$treatment=='CPS')

outliers_ache
dim(rvm_gsh)
rvm_gsh_ache <- rvm_gsh[-outliers_ache,]
dim(rvm_gsh_ache)
#View(rvm_gsh_ache)

###Welch Two Sample t-test
# Scenario 1
###24d ache with outliers
rvm_gsh$ache_ug_min_mg[which_controls]
#[1] 0.002543061 0.000725662 0.001163541 0.001621989 0.000737099 0.001184717 0.002179072 0.003789513 0.003437057 0.001573837
#[11] 0.002936909
rvm_gsh$ache_ug_min_mg[which_24d]
# [1] 0.001232294 0.001666101 0.002916909 0.004297287 0.001446049 0.001358612 0.004835887 0.002255579 0.001905794 0.002331805
# [11] 0.002349212
t.test(rvm_gsh$ache_ug_min_mg[which_controls],rvm_gsh$ache_ug_min_mg[which_24d])
#Welch Two Sample t-test
#data:  rvm_gsh$ache_ug_min_mg[which_controls] and rvm_gsh$ache_ug_min_mg[which_24d]
#t = -0.89072, df = 19.794, p-value = 0.3838
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.001429501  0.000574397
#sample estimates:
#  mean of x   mean of y 
#0.001990223 0.002417775 

# Scenario 2
###24d ache drop outliers
rvm_gsh_ache$ache_ug_min_mg[which_controls_ache]
#[1] 0.002543061 0.000725662 0.001163541 0.001621989 0.000737099 0.001184717 0.002179072 0.003789513 0.003437057 0.001573837
#[11] 0.002936909
rvm_gsh_ache$ache_ug_min_mg[which_24d_ache]
# [1] 0.001232294 0.001666101 0.002916909 0.001446049 0.001358612 0.002255579 0.001905794 0.002331805 0.002349212
t.test(rvm_gsh_ache$ache_ug_min_mg[which_controls_ache],rvm_gsh_ache$ache_ug_min_mg[which_24d_ache])
#data:  rvm_gsh_ache$ache_ug_min_mg[which_controls_ache] and rvm_gsh_ache$ache_ug_min_mg[which_24d_ache]
#t = 0.13418, df = 15.691, p-value = 0.895
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0007406310  0.0008405544
#sample estimates:
#  mean of x   mean of y 
#0.001990223 0.001940262 

### cohen d for 24d but drop outliers (no outliers)
# treatment levels first argument for this function
cohen.d(rvm_gsh_ache$ache_ug_min_mg[which_24d_ache], 
        rvm_gsh_ache$ache_ug_min_mg[which_controls_ache], 
        hedges.correction = T)
#Hedges's g
#g estimate: -0.05442714 (negligible)
#95 percent confidence interval:
#     lower      upper 
#-0.9589874  0.8501331 

###Welch Two Sample t-test logged
#with outliers logged
t.test(log(rvm_gsh$ache_ug_min_mg[which_controls]),log(rvm_gsh$ache_ug_min_mg[which_24d]))
#Welch Two Sample t-test

#data:  log(rvm_gsh$ache_ug_min_mg[which_controls]) and log(rvm_gsh$ache_ug_min_mg[which_24d])
#t = -1.099, df = 18.782, p-value = 0.2856
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.7051724  0.2198365
#sample estimates:
#  mean of x mean of y 
#-6.363308 -6.120640 




#######Chlorpyrifos
###Welch Two Sample t-test
###ache chlorpyrifos with outliers (no outliers so equivalent)
rvm_gsh$ache_ug_min_mg[which_controls]
#[1] 0.002543061 0.000725662 0.001163541 0.001621989 0.000737099 0.001184717 0.002179072 0.003789513 0.003437057 0.001573837
#[11] 0.002936909
rvm_gsh$ache_ug_min_mg[which_cps]
#[1] 0.003184610 0.001377392 0.002453795 0.002989280 0.002942805 0.001019211 0.002156165 0.003350097 0.001095101 0.001246339
#[11] 0.003166501
t.test(rvm_gsh$ache_ug_min_mg[which_controls],rvm_gsh$ache_ug_min_mg[which_cps])
#data:  rvm_gsh$ache_ug_min_mg[which_controls] and rvm_gsh$ache_ug_min_mg[which_cps]
#t = -0.65867, df = 19.624, p-value = 0.5178
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0011711762  0.0006095691
#sample estimates:
#  mean of x   mean of y 
#0.001990223 0.002271027 
### cohen d for cps but drop outliers (no outliers)
# treatment levels first argument for this function
cohen.d(rvm_gsh$ache_ug_min_mg[which_cps], 
        rvm_gsh$ache_ug_min_mg[which_controls], 
        hedges.correction = T)
#Hedges's g
#g estimate: 0.2701941 (small)
#95 percent confidence interval:
#     lower      upper 
#-0.5893823  1.1297704 

#drop outliers
t.test(rvm_gsh15$ache_ug_min_mg[which15_controls],rvm_gsh15$ache_ug_min_mg[which15_cps])


### LOGGED results
###Welch Two Sample t-test logged
###1:5 chlorpyrifos
#with or without outliers logged (outlier is 24d)
t.test(log(rvm_gsh$ache_ug_min_mg[which_controls]),log(rvm_gsh$ache_ug_min_mg[which_cps]))
#Welch Two Sample t-test

#data:  log(rvm_gsh$ache_ug_min_mg[which_controls]) and log(rvm_gsh$ache_ug_min_mg[which_cps])
#t = -0.81584, df = 19.191, p-value = 0.4246
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.6549764  0.2873999
#sample estimates:
#  mean of x mean of y 
#-6.363308 -6.179520 


#24d with outliers
t.test(log(rvm_gsh$ache_ug_min_mg[which_controls]),log(rvm_gsh$ache_ug_min_mg[which_24d]))
#Welch Two Sample t-test

#data:  log(rvm_gsh$ache_ug_min_mg[which_controls]) and log(rvm_gsh$ache_ug_min_mg[which_24d])
#t = -1.099, df = 18.782, p-value = 0.2856
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.7051724  0.2198365
#sample estimates:
#  mean of x mean of y 
#-6.363308 -6.120640 


###
#24d with logged outlier dropped
t.test(log(rvm_gsh$ache_ug_min_mg[which_logged_controls]),log(rvm_gsh$ache_ug_min_mg[which_logged_24d]))


