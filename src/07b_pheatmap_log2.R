dim(rvm_vst_mat)
summary(rvm_count_matrix)
rvm_count_matrix
pheatmap(log(rvm_count_matrix), annotation_col=rvm_metadata[,2:3], main="raw data")

colnames(rvm_count_matrix)
rvm_metadata[,2:3]

# sum of conncetrations by sample
rvm_sample_sums <- cbind(rvm_metadata[,2:3], t(rvm_count_matrix))
dim(rvm_sample_sums)
rvm_sample_sums$summed_concs <- rowSums(rvm_sample_sums[,3:29])
ggplot(rvm_sample_sums, aes(x=source, y=log(summed_concs), fill=treatment)) +
  geom_boxplot() +
  theme_bw() +
  labs(x="source", y="log(summed_concs)", fill="treatment") +
  theme(legend.position="top")


## normalization is a problem because we are assuming that the sum of these (significant metabolites)
## are equal across the samples
rvm_counts_normalized
dim(rvm_counts_normalized)
dim(rvm_metadata[,2:3])
pheatmap(log(rvm_counts_normalized), annotation_col=rvm_metadata[,2:3], main="normalized data")
colnames(rvm_counts_normalized)

pheatmap(log(rvm_counts_normalized[,c(1:11,34:44)]), annotation_col=rvm_metadata[,2:3], main="normalized data--Control")
pheatmap(log(rvm_counts_normalized[,c(12:22,45:55)]), annotation_col=rvm_metadata[,2:3], main="normalized data--Chlorpyrifos")
pheatmap(log(rvm_counts_normalized[,c(23:33,56:66)]), annotation_col=rvm_metadata[,2:3], main="normalized data--24D")

pheatmap(log(rvm_counts_normalized[,1:33]), annotation_col=rvm_metadata[,2:3], main="normalized data--Liver")
pheatmap(log(rvm_counts_normalized[,34:66]), annotation_col=rvm_metadata[,2:3], main="normalized data--Swabs")

dim(rvm_counts_normalized)
dim(rvm_metadata)

rvm_temp <- cbind(rvm_metadata, t(rvm_counts_normalized))
dim(rvm_temp)
colnames(rvm_temp)
rvm_means <- rvm_temp %>%
  group_by(treatment, source) %>%
  summarise(across(c(3:29), mean, na.rm=T))
rvm_means
combo_names <- paste0(rvm_means$source,"_",rvm_means$treatment)
dim(rvm_means)
dim(rvm_metadata)
View(rvm_means)
pheat_means <- t(log(rvm_means[,3:29]))
colnames(pheat_means) <- combo_names
pheat_means_annotate <- t(rvm_means[,1:2])
colnames(pheat_means_annotate) <- combo_names
pheatmap(pheat_means, main="means")
dim(t(log(rvm_means[,3:29])))
pheatmap(pheat_means, annotation_col=pheat_means_annotate, main="means")

pheat_vst_means <- t(log(rvm_means[,3:29]))
colnames(pheat_vst_means) <- combo_names
pheatmap(pheat_vst_means, main="transformed logged data (deseq)")


## liver
colnames(rvm_counts_normalized)
liver_cpf_fc1 <- rvm_counts_normalized[,12]/rvm_counts_normalized[,1]
liver_cpf_fc2 <- rvm_counts_normalized[,15]/rvm_counts_normalized[,4]
liver_cpf_fc3 <- rvm_counts_normalized[,16]/rvm_counts_normalized[,5]
liver_cpf_fc4 <- rvm_counts_normalized[,17]/rvm_counts_normalized[,6]
liver_cpf_fc5 <- rvm_counts_normalized[,18]/rvm_counts_normalized[,7]
liver_cpf_fc6 <- rvm_counts_normalized[,19]/rvm_counts_normalized[,8]
liver_cpf_fc7 <- rvm_counts_normalized[,20]/rvm_counts_normalized[,9]
liver_cpf_fc8 <- rvm_counts_normalized[,21]/rvm_counts_normalized[,10]
liver_cpf_fc9 <- rvm_counts_normalized[,22]/rvm_counts_normalized[,11]
liver_cpf_fc10 <- rvm_counts_normalized[,13]/rvm_counts_normalized[,2]
liver_cpf_fc11 <- rvm_counts_normalized[,14]/rvm_counts_normalized[,3]

liver_cpf_fc <- cbind(liver_cpf_fc1, liver_cpf_fc2, liver_cpf_fc3, liver_cpf_fc4, liver_cpf_fc5, liver_cpf_fc6,
      liver_cpf_fc7, liver_cpf_fc8, liver_cpf_fc9, liver_cpf_fc10, liver_cpf_fc11)
View(liver_cpf_fc)
pheatmap(liver_cpf_fc)

liver_24d_fc1 <- rvm_counts_normalized[,23]/rvm_counts_normalized[,1]
liver_24d_fc2 <- rvm_counts_normalized[,25]/rvm_counts_normalized[,4]
liver_24d_fc3 <- rvm_counts_normalized[,26]/rvm_counts_normalized[,5]
liver_24d_fc4 <- rvm_counts_normalized[,27]/rvm_counts_normalized[,6]
liver_24d_fc5 <- rvm_counts_normalized[,29]/rvm_counts_normalized[,7]
liver_24d_fc6 <- rvm_counts_normalized[,30]/rvm_counts_normalized[,8]
liver_24d_fc7 <- rvm_counts_normalized[,31]/rvm_counts_normalized[,9]
liver_24d_fc8 <- rvm_counts_normalized[,32]/rvm_counts_normalized[,10]
liver_24d_fc9 <- rvm_counts_normalized[,33]/rvm_counts_normalized[,11]
liver_24d_fc10 <- rvm_counts_normalized[,24]/rvm_counts_normalized[,2]
liver_24d_fc11 <- rvm_counts_normalized[,28]/rvm_counts_normalized[,3]

liver_24d_fc <- cbind(liver_24d_fc1, liver_24d_fc2, liver_24d_fc3, liver_24d_fc4, liver_24d_fc5, liver_24d_fc6,
                      liver_24d_fc7, liver_24d_fc8, liver_24d_fc9, liver_24d_fc10, liver_24d_fc11)
View(liver_24d_fc)
pheatmap(liver_24d_fc)

liver_fc <- cbind(liver_cpf_fc, liver_24d_fc)
pheatmap(log2(liver_fc))

##swab
colnames(rvm_counts_normalized)
swab_cpf_fc1 <- rvm_counts_normalized[,45]/rvm_counts_normalized[,34]
swab_cpf_fc2 <- rvm_counts_normalized[,48]/rvm_counts_normalized[,37]
swab_cpf_fc3 <- rvm_counts_normalized[,49]/rvm_counts_normalized[,38]
swab_cpf_fc4 <- rvm_counts_normalized[,50]/rvm_counts_normalized[,39]
swab_cpf_fc5 <- rvm_counts_normalized[,51]/rvm_counts_normalized[,40]
swab_cpf_fc6 <- rvm_counts_normalized[,52]/rvm_counts_normalized[,41]
swab_cpf_fc7 <- rvm_counts_normalized[,53]/rvm_counts_normalized[,42]
swab_cpf_fc8 <- rvm_counts_normalized[,54]/rvm_counts_normalized[,43]
swab_cpf_fc9 <- rvm_counts_normalized[,55]/rvm_counts_normalized[,44]
swab_cpf_fc10 <- rvm_counts_normalized[,46]/rvm_counts_normalized[,35]
swab_cpf_fc11 <- rvm_counts_normalized[,47]/rvm_counts_normalized[,36]

swab_cpf_fc <- cbind(swab_cpf_fc1, swab_cpf_fc2, swab_cpf_fc3, swab_cpf_fc4, swab_cpf_fc5, swab_cpf_fc6,
                      swab_cpf_fc7, swab_cpf_fc8, swab_cpf_fc9, swab_cpf_fc10, swab_cpf_fc11)
View(swab_cpf_fc)
pheatmap(log2(swab_cpf_fc))

swab_24d_fc1 <- rvm_counts_normalized[,56]/rvm_counts_normalized[,34]
swab_24d_fc2 <- rvm_counts_normalized[,59]/rvm_counts_normalized[,37]
swab_24d_fc3 <- rvm_counts_normalized[,66]/rvm_counts_normalized[,38]
swab_24d_fc4 <- rvm_counts_normalized[,60]/rvm_counts_normalized[,39]
swab_24d_fc5 <- rvm_counts_normalized[,61]/rvm_counts_normalized[,40]
swab_24d_fc6 <- rvm_counts_normalized[,62]/rvm_counts_normalized[,41]
swab_24d_fc7 <- rvm_counts_normalized[,63]/rvm_counts_normalized[,42]
swab_24d_fc8 <- rvm_counts_normalized[,64]/rvm_counts_normalized[,43]
swab_24d_fc9 <- rvm_counts_normalized[,65]/rvm_counts_normalized[,44]
swab_24d_fc10 <- rvm_counts_normalized[,57]/rvm_counts_normalized[,35]
swab_24d_fc11 <- rvm_counts_normalized[,58]/rvm_counts_normalized[,36]

swab_24d_fc <- cbind(swab_24d_fc1, swab_24d_fc2, swab_24d_fc3, swab_24d_fc4, swab_24d_fc5, swab_24d_fc6,
                      swab_24d_fc7, swab_24d_fc8, swab_24d_fc9, swab_24d_fc10, swab_24d_fc11)
View(swab_24d_fc)
pheatmap(log2(swab_24d_fc))

swab_fc <- cbind(swab_cpf_fc, swab_24d_fc)
pheatmap(log2(swab_fc))


fc_data <- data.frame(liver_fc, swab_fc)
fc_data
fc_tissue <- as.factor(c(rep("liver",22),rep("swab",22)))
fc_chemical <- as.factor(c(rep("cpf",11),rep("24d",11),rep("cpf",11),rep("24d",11)))
fc_annotate <- data.frame(fc_tissue, fc_chemical)
fc_annotate
dim(fc_annotate)
rownames(fc_annotate) <- colnames(fc_data)
dim(fc_data)
pheatmap(log2(fc_data), main="fold change")
pheatmap(log2(fc_data), annotation_col=fc_annotate, main="fold change")

# means of fold changes
swab_24d_fc
liver_cpf_fc_means <- rowMeans(liver_cpf_fc)
liver_24d_fc_means <- rowMeans(liver_24d_fc)
swab_cpf_fc_means <- rowMeans(swab_cpf_fc)
swab_24d_fc_means <- rowMeans(swab_24d_fc)

fc_data_means <- cbind(liver_cpf_fc_means, liver_24d_fc_means, swab_cpf_fc_means, swab_24d_fc_means)
pheatmap(log2(fc_data_means), main="fold change means")
pheatmap(log2(fc_data_means[-15,]), main="fold change means")
