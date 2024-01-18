# Van Meter salamander stat analysis
# 24d and chlorpyrifos as stressors/exposures
# responses are glutathione and acetylcholinesterase


#Install and load supporting libraries.
print(Sys.info()[4])
#BiocManager::install("pathview")

library(dplyr)
library(forcats)
library(FSA)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(BGGM)
library(ggm)
library(corrplot)
library(qgraph)
library(robFitConGraph)
library(reshape2)
library(glasso)
library(igraph)
library(GGally)
library(matrixcalc)
library(Matrix)
library(pathview)
library(pheatmap)
library(effsize)

#browseVignettes("pathview")

print("list of loaded packages: ")
print((.packages()))

#tom epa windows
if(Sys.info()[4]=="DZ2626UTPURUCKE"){
  rvm_root <- file.path("c:", "git", "vanmeter_salamander_exposures")
}
#tom epa windows 2
if(Sys.info()[4]=="LZ26TPURUCKE-2"){
  rvm_root <- file.path("c:", "Users", "tpurucke", "git", "vanmeter_salamander_exposures")
}

print(paste("Root directory location: ", rvm_root, sep=""))

rvm_data_in <- file.path(rvm_root, "data_in")
rvm_data_out <- file.path(rvm_root, "data_out")
rvm_graphics <- file.path(rvm_root, "graphics")

#check to see if directories are accessible
boo = file.exists(file.path(rvm_data_in,"Final_GSH_Salamanders_April_2020.csv"))
print(paste("check to see if R can access GSF file OK: ", boo))

#GSH units = nM/mL
#Weight (g)
#SVL (mm)

#gsh data (ache added 5/2022) # liver data
rvm_gsh <- read.csv(file.path(rvm_data_in,"/Final_GSH_Salamanders_April_2020.csv"), stringsAsFactors = TRUE)
#View(rvm_gsh)
colnames(rvm_gsh)[1] <- 'treatment'

#gsh data swabs
rvm_gsh_swabs <- read.csv(file.path(rvm_data_in,"/gsh_swab_concs.csv"), stringsAsFactors = TRUE)
#View(rvm_gsh_swabs)

# Using the average of gsh 1_5 and 1_8 for manuscript
rvm_gsh$gsh_nM_mL <- rowMeans(cbind(rvm_gsh$gsh_1_5_dilution_nM_mL, rvm_gsh$gsh_1_8_dilution_nM_mL))

summary(rvm_gsh)
colnames(rvm_gsh)
#View(rvm_gsh)

# peaks from wmh/dag
# metabolites with multiple peas were summed into one peak

# combined_enrichment
rvm_combined_enrichment <- read.csv(file.path(rvm_data_in,"/combined_enrichment.csv"), stringsAsFactors = TRUE)
colnames(rvm_combined_enrichment)
rownames(rvm_combined_enrichment)
dim(rvm_combined_enrichment)

# all peaks
rvm_all_peaks <- read.csv(file.path(rvm_data_in,"/all_significant_peaks.csv"), stringsAsFactors = TRUE)
colnames(rvm_all_peaks)
rownames(rvm_all_peaks)
dim(rvm_all_peaks)

# urea cycle
rvm_urea <- read.csv(file.path(rvm_data_in,"/urea_cycle_peaks.csv"), stringsAsFactors = TRUE)
colnames(rvm_urea)
rownames(rvm_urea)
dim(rvm_urea)

# gly ser metabolism
rvm_gly_ser_metabolism <- read.csv(file.path(rvm_data_in,"/gly_ser_metabolism_peaks.csv"), stringsAsFactors = TRUE)
dim(rvm_gly_ser_metabolism)

# glu metabolism
rvm_glu_metabolism <- read.csv(file.path(rvm_data_in,"/glu_metabolism_peaks.csv"), stringsAsFactors = TRUE)
dim(rvm_glu_metabolism)

# gluconeogenesis
rvm_gluconeogenesis <- read.csv(file.path(rvm_data_in,"/gluconeogenesis_peaks.csv"), stringsAsFactors = TRUE)
dim(rvm_gluconeogenesis)

# glu ala cycle
rvm_glu_ala_cycle <- read.csv(file.path(rvm_data_in,"/glu_ala_cycle_peaks.csv"), stringsAsFactors = TRUE)
dim(rvm_glu_ala_cycle)

# glu ala cycle
cpf_ratios <- read.csv(file.path(rvm_data_in,"cpf_ratios.csv"), stringsAsFactors = TRUE)
dim(cpf_ratios)


