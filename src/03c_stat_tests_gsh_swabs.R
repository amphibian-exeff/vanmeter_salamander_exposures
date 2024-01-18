############
# not an anova since no interaction, do ttests
# indices for the 3 scenarios: all data, dropping outlier, logged with dropped outlier
# all data
which_swab_controls <- which(rvm_gsh_swabs$treatment=='CON')
which_swab_24d <- which(rvm_gsh_swabs$treatment=='D')
which_swab_cps <- which(rvm_gsh_swabs$treatment=='CHL')
# dropping outliers for averaged dilutions
dim(rvm_gsh_swabs_drop_outliers)
which_dropout_swab_controls <- which(rvm_gsh_swabs_drop_outliers$treatment=='CON')
which_dropout_swab_24d <- which(rvm_gsh_swabs_drop_outliers$treatment=='D')
which_dropout_swab_cps <- which(rvm_gsh_swabs_drop_outliers$treatment=='CHL')



####### 24D SWABS
### these are used in the manuscript for swabs
###
### gsh swab t-test for 24d but drop outliers
rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls]
#  [1] 0.1460550 0.1351652 0.1282957 0.1169261 0.1296864 0.1386652 0.1297304 0.1443000 0.1186521 0.1380652
rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_24d]
#  [1] 0.1683391 0.1495550 0.1371304 0.1211692 0.1627043 0.1215609 0.1241609 0.1065349 0.1453000 0.1472000
t.test(rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls], 
       rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_24d])
#data:  rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls] and rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_24d]
#t = -0.83082, df = 13.161, p-value = 0.4209
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.020903766  0.009281065
#sample estimates:
#  mean of x mean of y 
#0.1325541 0.1383655 

### cohen d for 24d but drop outliers
# treatment levels first argument for this function
cohen.d(rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_24d], 
        rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls], 
        hedges.correction = T)
#g estimate: 0.3558548 (small)
#95 percent confidence interval:
#  lower      upper 
#-0.5511004  1.2628100 


### gsh swab t-test for chlorpyrifos but drop outliers
rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls]
#  [1] 0.1460550 0.1351652 0.1282957 0.1169261 0.1296864 0.1386652 0.1297304 0.1443000 0.1186521 0.1380652
rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_cps]
#  0.1569696 0.1567696 0.1238692 0.1339036 0.1326957 0.1319957 0.1789237 0.1461000 0.1276609
t.test(rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls], 
       rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_cps])
#data:  rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls] and rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_cps]
#t = -1.5766, df = 12.078, p-value = 0.1407
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.025371213  0.004059975
#sample estimates:
#  mean of x mean of y 
#0.1325541 0.1432098 

### cohens d for chlorpyrifos but drop outliers
cohen.d(rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_cps],
        rvm_gsh_swabs_drop_outliers$total_GSH[which_dropout_swab_controls], 
        hedges.correction = T)
# Hedges's g
#g estimate: 0.7132158 (medium)
#95 percent confidence interval:
#     lower      upper 
#-0.2416784  1.6681101  
