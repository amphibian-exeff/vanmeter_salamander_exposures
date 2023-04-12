dim(rvm_gsh)
summary(rvm_gsh)
colnames(rvm_gsh)

#treatments
unique(rvm_gsh$treatment)
#[1] C   24D CPS
#Levels: 24D C CPS

## summary stats, normality tests, and check for outliers for gsh



####################
Summarize(gsh_1_5_dilution_nM_mL ~ treatment, data=rvm_gsh, digits=3)
#treatment  n   mean    sd   min     Q1 median    Q3    max
#1       24D 11 14.696 5.347 4.693 11.467 15.461 17.67 24.971
#2         C 11 12.239 6.056 5.962  8.139 11.640 14.14 27.479
#3       CPS 11 10.460 4.003 6.862  7.657  8.232 12.85 18.832

Summarize(gsh_1_8_dilution_nM_mL ~ treatment, data=rvm_gsh, digits=3)
#treatment  n   mean    sd   min     Q1 median     Q3    max
#1       24D 11 17.982 5.352 5.831 16.028 19.246 21.693 23.731
#2         C 11 15.335 7.881 7.898 10.367 14.663 17.284 36.158
#3       CPS 11 13.610 4.608 7.135 11.088 12.337 15.818 24.071

Summarize(gsh_nM_mL ~ treatment, data=rvm_gsh, digits=3)
#treatment  n   mean    sd   min     Q1 median     Q3    max
#1       24D 11 16.339 5.133 5.262 13.748 17.354 20.495 22.655
#2         C 11 13.787 6.949 6.954  9.253 13.045 15.712 31.818
#3       CPS 11 12.035 4.255 6.999  9.520 10.255 14.334 21.452

# normality tests
# shapiro.test
# p.value: an approximate p-value for the test. This is said in Royston (1995) 
# to be adequate for p.value < 0.1.
aggregate(gsh_nM_mL ~ treatment,
          data = rvm_gsh,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment gsh_nM_mL.W gsh_nM_mL.V2
#1       24D  0.93062706   0.41723177
#2         C  0.81890542   0.01668438
#3       CPS  0.89189592   0.14684800

aggregate(gsh_1_5_dilution_nM_mL ~ treatment,
          data = rvm_gsh,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment gsh_1_5_dilution_nM_mL.W gsh_1_5_dilution_nM_mL.V2
#1       24D               0.97338920                0.91853993
#2         C               0.84681669                0.03878205
#3       CPS               0.84809315                0.04030010

aggregate(gsh_1_8_dilution_nM_mL ~ treatment,
          data = rvm_gsh,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment gsh_1_8_dilution_nM_mL.W gsh_1_8_dilution_nM_mL.V2
#1       24D              0.901545590               0.193046103
#2         C              0.793815305               0.007799802
#3       CPS              0.918494639               0.306301002

#identify outliers
colnames(rvm_gsh)
head(rvm_gsh)
# identify outliers
rvm_gsh %>%
  group_by(treatment) %>%
  identify_outliers("gsh_nM_mL")
#treatment salamander_id weight_g svl_mm m_f   gsh_1_5_dilution… gsh_1_8_dilutio… ache_ug_min_mg log_gsh_15 log_gsh_18
#<fct>     <fct>            <dbl>  <dbl> <fct>             <dbl>            <dbl>          <dbl>      <dbl>      <dbl>
#  1 C         CONS3             22.0   93.3 M                  27.5             36.2        0.00116       3.31       3.59

rvm_gsh %>%
  group_by(treatment) %>%
  identify_outliers("gsh_1_5_dilution_nM_mL")
#treatment salamander_id weight_g svl_mm m_f   gsh_1_5_dilution_nM_mL gsh_1_8_dilution_nM_mL is.outlier is.extreme
#<fct>     <fct>            <dbl>  <dbl> <fct>                  <dbl>                  <dbl> <lgl>      <lgl>     
#  1 C         CONS3             22.0   93.3 M                       27.5                   36.2 TRUE       FALSE

rvm_gsh %>%
  group_by(treatment) %>%
  identify_outliers("gsh_1_8_dilution_nM_mL")
#treatment salamander_id weight_g svl_mm m_f   gsh_1_5_dilution_nM_mL gsh_1_8_dilution_nM_mL is.outlier is.extreme
#<fct>     <fct>            <dbl>  <dbl> <fct>                  <dbl>                  <dbl> <lgl>      <lgl>     
#  1 24D       DS7               17.8   87.2 M                       4.69                   5.83 TRUE       FALSE     
#2 C         CONS3             22.0   93.3 M                      27.5                   36.2  TRUE       FALSE     
#3 CPS       CHLS4             15.3   82.2 M                      18.8                   24.1  TRUE       FALSE   

outliers_158 <- which(rvm_gsh$salamander_id=='CONS3')
outliers_158
outliers15 <- which(rvm_gsh$salamander_id=='CONS3')
outliers15
outliers18 <- c(which(rvm_gsh$salamander_id=='CONS3'),which(rvm_gsh$salamander_id=='DS7'),
                which(rvm_gsh$salamander_id=='CHLS4'))
outliers18

# repeat normality tests but drop outliers
dim(rvm_gsh)
rvm_gsh_158 <- rvm_gsh[-outliers_158,]
dim(rvm_gsh_158)
aggregate(gsh_nM_mL ~ treatment,
          data = rvm_gsh_158,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment gsh_nM_mL.W gsh_nM_mL.V2
#1       24D   0.9306271    0.4172318
#2         C   0.9394081    0.5464464
#3       CPS   0.8918959    0.1468480

rvm_gsh15 <- rvm_gsh[-outliers15,]
dim(rvm_gsh_158)
aggregate(gsh_1_5_dilution_nM_mL ~ treatment,
          data = rvm_gsh15,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment gsh_1_5_dilution_nM_mL.W gsh_1_5_dilution_nM_mL.V2
#1       24D                0.9733892                 0.9185399
#2         C                0.9297638                 0.4455725
#3       CPS                0.8480932                 0.0403001

rvm_gsh18 <- rvm_gsh[-outliers18,]
dim(rvm_gsh18)
aggregate(gsh_1_8_dilution_nM_mL ~ treatment,
          data = rvm_gsh18,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment gsh_1_8_dilution_nM_mL.W gsh_1_8_dilution_nM_mL.V2
#1       24D                0.9480377                 0.6453436
#2         C                0.9303179                 0.4510260
#3       CPS                0.9586766                 0.7706828

##### logged
### normality tests on logged data
aggregate(log(gsh_nM_mL) ~ treatment,
          data = rvm_gsh,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment log(gsh_nM_mL).W log(gsh_nM_mL).V2
#1       24D       0.81116854        0.01319703
#2         C       0.94720449        0.60846772
#3       CPS       0.95392530        0.69435838

aggregate(log(gsh_1_5_dilution_nM_mL) ~ treatment,
          data = rvm_gsh,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment log(gsh_1_5_dilution_nM_mL).W log(gsh_1_5_dilution_nM_mL).V2
#1       24D                     0.8915908                      0.1455693
#2         C                     0.9476519                      0.6141094
#3       CPS                     0.8873104                      0.1286907

aggregate(log(gsh_1_8_dilution_nM_mL) ~ treatment,
          data = rvm_gsh,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment log(gsh_1_8_dilution_nM_mL).W log(gsh_1_8_dilution_nM_mL).V2
#1       24D                   0.771880265                    0.004017602
#2         C                   0.939156563                    0.510633852
#3       CPS                   0.971349290                    0.899813724

#add logged columns
rvm_gsh$log_gsh_158 <- log(rvm_gsh$gsh_nM_mL)
rvm_gsh$log_gsh_15 <- log(rvm_gsh$gsh_1_5_dilution_nM_mL)
rvm_gsh$log_gsh_18 <- log(rvm_gsh$gsh_1_8_dilution_nM_mL)

#identify outliers in logged data
# identify outliers
rvm_gsh %>%
  group_by(treatment) %>%
  identify_outliers("log_gsh_158")
#treatment salamander_id weight_g svl_mm m_f   gsh_1_5_dilution… gsh_1_8_dilutio… ache_ug_min_mg log_gsh_15 log_gsh_18
#<fct>     <fct>            <dbl>  <dbl> <fct>             <dbl>            <dbl>          <dbl>      <dbl>      <dbl>
#  1 24D       DS7               17.8   87.2 M                  4.69             5.83        0.00484       1.55       1.76

rvm_gsh %>%
  group_by(treatment) %>%
  identify_outliers("log_gsh_15")
# A tibble: 1 x 11
#treatment salamander_id weight_g svl_mm m_f   gsh_1_5_dilutio~ gsh_1_8_dilutio~ log_gsh_15 log_gsh_18 is.outlier
#<fct>     <fct>            <dbl>  <dbl> <fct>            <dbl>            <dbl>      <dbl>      <dbl> <lgl>     
#  1 24D       DS7               17.8   87.2 M                 4.69             5.83       1.55       1.76 TRUE 
# identify outliers
rvm_gsh %>%
  group_by(treatment) %>%
  identify_outliers("log_gsh_18")
# A tibble: 1 x 11
#treatment salamander_id weight_g svl_mm m_f   gsh_1_5_dilutio~ gsh_1_8_dilutio~ log_gsh_15 log_gsh_18 is.outlier
#<fct>     <fct>            <dbl>  <dbl> <fct>            <dbl>            <dbl>      <dbl>      <dbl> <lgl>     
#  1 24D       DS7               17.8   87.2 M                 4.69             5.83       1.55       1.76 TRUE  

outliers_logged <- which(rvm_gsh$salamander_id=='DS7')
# repeat normality tests for logged data but drop outliers
rvm_gsh_logged <- rvm_gsh[-outliers_logged,]
### normality tests on logged data
aggregate(log(gsh_nM_mL) ~ treatment,
          data = rvm_gsh_logged,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment log(gsh_nM_mL).W log(gsh_nM_mL).V2
#1       24D        0.9344307         0.4928570
#2         C        0.9472045         0.6084677
#3       CPS        0.9539253         0.6943584

aggregate(log(gsh_1_5_dilution_nM_mL) ~ treatment,
          data = rvm_gsh_logged,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment log(gsh_1_5_dilution_nM_mL).W log(gsh_1_5_dilution_nM_mL).V2
#1       24D                     0.9548845                      0.7263417
#2         C                     0.9476519                      0.6141094
#3       CPS                     0.8873104                      0.1286907

aggregate(log(gsh_1_8_dilution_nM_mL) ~ treatment,
          data = rvm_gsh_logged,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment log(gsh_1_8_dilution_nM_mL).W log(gsh_1_8_dilution_nM_mL).V2
#1       24D                     0.9244898                      0.3959265
#2         C                     0.9391566                      0.5106339
#3       CPS                     0.9713493                      0.8998137



### SWAB data for glutathione
colnames(rvm_gsh_swabs)
Summarize(total_GSH ~ treatment, data=rvm_gsh_swabs, digits=3)
#treatment  n  mean    sd   min    Q1 median    Q3   max
#1       CHL 11 0.150 0.047 0.085 0.130  0.134 0.157 0.273
#2       CON 11 0.129 0.016 0.091 0.123  0.130 0.138 0.146
#3         D 11 0.148 0.038 0.107 0.123  0.145 0.156 0.247

# normality tests--swabs
# shapiro.test
# p.value: an approximate p-value for the test. This is said in Royston (1995) 
# to be adequate for p.value < 0.1.
aggregate(total_GSH ~ treatment,
          data = rvm_gsh_swabs,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment total_GSH.W total_GSH.V2
#1       CHL  0.81027986   0.01284627
#2       CON  0.87454060   0.08864277
#3         D  0.81657098   0.01554492

# identify outliers
rvm_gsh_swabs %>%
  group_by(treatment) %>%
  identify_outliers("total_GSH")
#treatment ID    i_slope y_intercept f_slope total_GSH is.outlier is.extreme
#<fct>     <fct>   <dbl>       <dbl>   <dbl>     <dbl> <lgl>      <lgl>     
#1 CHL       CHL6   0.0069       0.102  0.0184    0.273  TRUE       TRUE      
#2 CHL       CHL8   0.0032       0.105  0.0169    0.0847 TRUE       FALSE     
#3 CON       CON8   0.0033       0.105  0.0169    0.0907 TRUE       FALSE     
#4 D         D10    0.0062       0.120  0.0169    0.247  TRUE       FALSE  

# save swab outliers
outliers_swabs <- c(which(rvm_gsh_swabs$ID=='CHL6'), which(rvm_gsh_swabs$ID=='CHL8'),
                which(rvm_gsh_swabs$ID=='CON8'), which(rvm_gsh_swabs$ID=='D10'))
outliers_swabs

# drop outliers and rerun tests
# repeat normality tests but drop outliers
dim(rvm_gsh_swabs)
colnames(rvm_gsh_swabs_drop_outliers)
rvm_gsh_swabs_drop_outliers <- rvm_gsh_swabs[-outliers_swabs,]
dim(rvm_gsh_swabs_drop_outliers)
aggregate(total_GSH ~ treatment,
          data = rvm_gsh_swabs_drop_outliers,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#treatment total_GSH.W total_GSH.V2
#1       CHL   0.8911902    0.2052204
#2       CON   0.9450850    0.6108477
#3         D   0.9603279    0.7896053

