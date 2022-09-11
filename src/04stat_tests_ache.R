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
#with outliers
rvm_gsh$ache_ug_min_mg[which_controls]
rvm_gsh$ache_ug_min_mg[which_24d]
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
###1:5 chlorpyrifos
#with outliers
t.test(rvm_gsh$ache_ug_min_mg[which_controls],rvm_gsh$ache_ug_min_mg[which_cps])
#Welch Two Sample t-test

#data:  rvm_gsh$ache_ug_min_mg[which_controls] and rvm_gsh$ache_ug_min_mg[which_cps]
#t = -0.65867, df = 19.624, p-value = 0.5178
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0011711762  0.0006095691
#sample estimates:
#  mean of x   mean of y 
#0.001990223 0.002271027 


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


