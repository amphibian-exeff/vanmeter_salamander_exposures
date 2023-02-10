dim(cpf_ratios)
colnames(cpf_ratios)
rownames(cpf_ratios)

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
cpf_ratios_log2 <- log2(cpf_ratios[,3:13])
colnames(cpf_ratios_log2)
rownames(cpf_ratios_log2)<- kegg_list
cpf_ratios_log2
# kegg pathway ids
pathview(cpd.data=cpf_ratios_log2, pathway.id="00640") # pathway id 5 numbers "00640"
# Xenopus laevis (xla), tropicalis (xtr)
#pathview(cpd.data=cpf_ratios_log2, pathway.id="04612", pathway.name="hsa04612") # pathway name 3 letter sp then numbers ""

# urea cycle = 00220
pathview(cpd.data=cpf_ratios_log2, pathway.id="00220")
# glycine serine metabolism
# glucose metabolism
# gluconeogenesis
# glucose alanine cycle
