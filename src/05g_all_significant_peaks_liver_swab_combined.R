
# heatmap of normalized metabolite concentrations that include the controls
pheatmap(rvm_all_peaks_heatmap_values2_liver)
pheatmap(rvm_all_peaks_heatmap_values2_swab)

colnames(rvm_all_peaks_heatmap_values2_liver)
colnames(rvm_all_peaks_heatmap_values2_swab)

rvm_all_peaks_heatmap_values2_combined <- cbind(rvm_all_peaks_heatmap_values2_liver,rvm_all_peaks_heatmap_values2_swab)
colnames(rvm_all_peaks_heatmap_values2_combined) <- c("x24D_liver", "Chlorpyrifos_liver", "Control_liver", 
                                                      "x24D_swab", "Chlorpyrifos_swab", "Control_swab" )

pheatmap(rvm_all_peaks_heatmap_values2_combined)


# heatmap of log2 fold change between treatments and respective controls
# original concentration matrices
dim(rvm_all_peaks)
colnames(rvm_all_peaks)
rownames(rvm_all_peaks)
rvm_all_peaks$sample
rvm_all_peaks$treatment

df_rvm_peaks_con <- rvm_all_peaks %>%
  filter(treatment=="Control")
df_rvm_peaks_con$sample

df_rvm_peaks_24D <- rvm_all_peaks %>%
  filter(treatment=="24D")
df_rvm_peaks_con$sample

df_rvm_peaks_chlrpyr <- rvm_all_peaks %>%
  filter(treatment=="Chlorpyrifos")
