#FDR values

dim(rvm_combined_enrichment)
summary(rvm_combined_enrichment)
rownames(rvm_combined_enrichment) <- rvm_combined_enrichment$metabolite_set
rvm_combined_enrichment2 <- rvm_combined_enrichment[,-1]
pheatmap(-log(rvm_combined_enrichment2))
