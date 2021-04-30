#Install and load supporting libraries.
print(Sys.info()[4])

library(dplyr)
library(forcats)
library(FSA)
library(ggplot2)
library(ggpubr)
library(rstatix)

print("list of loaded packages: ")
print((.packages()))

#tom epa windows
if(Sys.info()[4]=="DZ2626UTPURUCKE"){
  rvm_root <- file.path("c:", "git", "rvm_salamander_gsh")
}
if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  rvm_root <- file.path("c:","git","rvm_salamander_gsh")
}

print(paste("Root directory location: ", rvm_root, sep=""))

rvm_data_in <- file.path(rvm_root, "data_in")
rvm_data_out <- file.path(rvm_root, "data_out")
rvm_graphics <- file.path(rvm_root, "graphics")

#check to see if directories are accessible
boo = file.exists(file.path(rvm_data_in,"/Final_GSH_Salamanders_April_2020.csv"))
print(paste("check to see if R can access GSF file OK: ", boo))

#GSH units = nM/mL
#Weight (g)
#SVL (mm)

rvm_gsh <- read.csv(file.path(rvm_csv_in,"/Final_GSH_Salamanders_April_2020.csv"), stringsAsFactors = TRUE)
View(rvm_gsh)
colnames(rvm_gsh)[1] <- 'treatment'

summary(rvm_gsh)
colnames(rvm_gsh)
View(rvm_gsh)

