dim(cpf_ratios)
colnames(cpf_ratios) # columns are individuals, plus kegg and common metabolite name
rownames(cpf_ratios) # rows are metabolites

# switch from cpf_ratios to straigth averages before finding the ratios

metabolite_list <- cpf_ratios$sample
metabolite_list
kegg_list <- cpf_ratios$kegg_compound
kegg_list

cpf_ratios_t <- t(log2(cpf_ratios[,3:13]))
rownames(cpf_ratios_t)
colnames(cpf_ratios_t) <- metabolite_list
as.numeric(cpf_ratios_t)
cpf_ratios_t <- as.data.frame(cpf_ratios_t)

library(purrr)
library(tidyr)
library(ggplot2)

cpf_ratios_t %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

# for pathview
cpf_ratios_log2 <- as.data.frame(log2(rowMeans(cpf_ratios[,3:13])))
rowMeans(cpf_ratios[,3:13])
colnames(cpf_ratios_log2)
rownames(cpf_ratios_log2)<- kegg_list
cpf_ratios_log2

# kegg pathway ids
# proponoate metabolism = 00640 (gluconeogenesis)
#Propionate metabolism is the set of biochemical processes that break down the short-chain fatty acid propionate. 
#Propionate is a three-carbon fatty acid that is produced by gut bacteria during the fermentation of dietary fiber 
#in the digestive system, and it is also produced as a metabolic byproduct of certain amino acids and other metabolic processes.
#In humans and other animals, propionate is metabolized primarily in the liver via beta-oxidation. In beta-oxidation, 
#propionate is broken down into acetyl-CoA, which can then be used by the body as a source of energy or as a 
#precursor molecule for the synthesis of other compounds.
#One of the main products of propionate metabolism in the liver is glucose. Propionate is converted to succinate, 
#which is then converted to oxaloacetate and finally to glucose via gluconeogenesis. 
#his process is important for maintaining blood glucose levels and providing energy to the body's tissues.

pathview(cpd.data=cpf_ratios_log2, species='hsa', pathway.id="00640") # pathway id 5 numbers "00640"
pathview(cpd.data=cpf_ratios_log2, species='xla', pathway.id="00640")

# urea cycle = 00220
#Arginine biosynthesis is the process by which living organisms synthesize the amino acid arginine from 
#precursor molecules. Arginine is an essential amino acid in humans and many other organisms, obtained through the 
#diet or synthesized through biosynthesis pathways.
#In animals, including humans, arginine biosynthesis occurs through the citrulline-NO cycle. In this pathway, 
#citrulline is synthesized from ornithine by the enzyme ornithine transcarbamylase, and then citrulline is 
#converted to arginine by the enzyme argininosuccinate synthase. Arginine is then used in a variety of 
#metabolic processes throughout the body, including the synthesis of nitric oxide, creatine, and urea.
pathview(cpd.data=cpf_ratios_log2, species='hsa', pathway.id="00220")
pathview(cpd.data=cpf_ratios_log2, species='xla', pathway.id="00220")

# glycine serine threonine metabolism = 00260
pathview(cpd.data=cpf_ratios_log2, species='hsa', pathway.id="00260")
pathview(cpd.data=cpf_ratios_log2, species='xla', pathway.id="00260")

# gluconeogenesis, glucose metabolism = 00010
pathview(cpd.data=cpf_ratios_log2, species='hsa', pathway.id="00010")
pathview(cpd.data=cpf_ratios_log2, species='xla', pathway.id="00010")

# glucose alanine cycle = 00910
pathview(cpd.data=cpf_ratios_log2, species='hsa', pathway.id="00910")
pathview(cpd.data=cpf_ratios_log2, species='xla', pathway.id="00910")

