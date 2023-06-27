dim(rvm_all_peaks)
colnames(rvm_all_peaks)
unique(rvm_all_peaks$source)
unique(rvm_all_peaks$treatment)
unique(rvm_all_peaks$sample)

# metadata file, samples as rows
rvm_metadata <- rvm_all_peaks[,1:3]
rvm_metadata$sample2 <- paste(rvm_metadata$source,"_",rvm_metadata$sample, sep="")
rvm_metadata$sample2[28] <-"liver_d4_rep"
rownames(rvm_metadata) <- rvm_metadata$sample2
rownames(rvm_metadata)
unique(rvm_metadata$treatment)

# set reference level
#converting counts to integer mode
#it appears that the last variable in the design formula, 'treatment',
#has a factor level, 'Control', which is not the reference level. we recommend
#to use factor(...,levels=...) or relevel() to set this as the reference level
#before proceeding. for more information, please see the 'Note on factor levels'
#in vignette('DESeq2')
levels(rvm_metadata$treatment)

# count matrix, samples as columns
rvm_count_matrix <- t(round(rvm_all_peaks[,4:30]))
dim(rvm_count_matrix)
colnames(rvm_count_matrix) <- rvm_metadata$sample2
rownames(rvm_count_matrix)

# check to confirm equal
all(rownames(rvm_metadata)==colnames(rvm_count_matrix))

# create deseq object
library(DESeq2)
rvm_dsq2 <- DESeqDataSetFromMatrix(countData = rvm_count_matrix,
                                   colData=rvm_metadata,
                                   design = ~ source + treatment) # treatment is of interest so last

rvm_dsq2_interaction <- DESeqDataSetFromMatrix(countData = rvm_count_matrix,
                                   colData=rvm_metadata,
                                   design = ~ source + treatment + source:treatment) # adds the interaction term "effect of source on the effect of treatment"


# normalize raw counts
# fills a slot in the deseq2 object
rvm_dsq2 <- estimateSizeFactors(rvm_dsq2)
# divides by these size factors when normalizing
sizeFactors(rvm_dsq2)
rvm_counts_normalized <- counts(rvm_dsq2, normalized=T)
#View(rvm_counts_normalized)

# unsupervised clustering analysis
# log transformation via variance stabilizing transformation
# vst(rvm_counts_normalized) #error because too few metabolites
rvm_vst <- varianceStabilizingTransformation(rvm_dsq2, blind=T)
rvm_vst_mat <- assay(rvm_vst)
rvm_vst_cor <- cor(rvm_vst_mat)
rvm_vst_cor2 <- cor(rvm_counts_normalized)
rvm_vst_cor3 <- cor(rvm_count_matrix)
#View(rvm_vst_cor)
colnames(rvm_metadata)
pheatmap(rvm_vst_cor, annotation_col=rvm_metadata[,2:3], main="variance stabilizing transformation on logged data")
pheatmap(rvm_vst_cor2, annotation_col=rvm_metadata[,2:3], main="normalized count data")
pheatmap(rvm_vst_cor3, annotation_col=rvm_metadata[,2:3], main="raw data")

# principal components analysis
plotPCA(rvm_vst, intgroup=c("source","treatment"))

# run deseq--now contains everything
rvm_contrast <- DESeq(rvm_dsq2)
plotDispEsts(rvm_contrast)

# deseq results
rvm_swab_liver <- results(rvm_contrast, 
        contrast = c("source", "swab", "liver"),
        alpha = 0.05)
plotMA(rvm_swab_liver, ylim=c(-5,5))
mcols(rvm_swab_liver)
# Generate logical column 
smoc2_res_all <- data.frame(smoc2_res) %>% mutate(threshold = padj < 0.05)