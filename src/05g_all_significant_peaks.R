#dim(rvm_all_peaks)
dim(rvm_all_peaks)
rvm_all_peaks <- rvm_all_peaks
#View(rvm_all_peaks)

rvm_all_peaks_mat <- rvm_all_peaks[,3:ncol(rvm_all_peaks)]
dim(rvm_all_peaks_mat)
ncols_all_peaks <- dim(rvm_all_peaks_mat)[[2]]
#View(rvm_all_peaks_mat)
rvm_all_peaks_mat_log <- log2(rvm_all_peaks_mat)

#View(rvm_all_peaks_mat_log)
# standardize is in 00support_functions.R
rvm_all_peaks_std2 = normal_score_transform(rvm_all_peaks_mat, 1:ncols_all_peaks)
rvm_all_peaks_std2 <- cbind(rvm_all_peaks[,2],rvm_all_peaks_std)
colnames(rvm_all_peaks_std2)[[1]] <- "treatment"
#View(rvm_all_peaks_std) # standardized after log2 transformation

# summarize for a heatmap
dim(rvm_all_peaks_std2)
colnames(rvm_all_peaks_std2)
rvm_all_peaks_heatmap_values <- rvm_all_peaks_std2 %>%
  group_by(treatment) %>%
  summarise_all(.funs = c(mean="mean"))
  #summarise_all(across(starts_with("numeric_column"), mean, na.rm = TRUE))

rvm_all_peaks_heatmap_values2 <- t(rvm_all_peaks_heatmap_values)
colnames(rvm_all_peaks_heatmap_values2) <- c("x24D","Chlorpyrifos","Control")
rvm_all_peaks_heatmap_values2 <- as.data.frame(rvm_all_peaks_heatmap_values2[-1,])
rvm_all_peaks_heatmap_values2$x24D <- as.numeric(rvm_all_peaks_heatmap_values2$x24D)
rvm_all_peaks_heatmap_values2$Chlorpyrifos <- as.numeric(rvm_all_peaks_heatmap_values2$Chlorpyrifos)
rvm_all_peaks_heatmap_values2$Control <- as.numeric(rvm_all_peaks_heatmap_values2$Control)
rownames(rvm_all_peaks_heatmap_values2) <- gsub("_mean", "", rownames(rvm_all_peaks_heatmap_values2))
pheatmap(rvm_all_peaks_heatmap_values2)

rvm_all_peaks_std_control <- rvm_all_peaks_std[1:11,]
#View(rvm_all_peaks_std_control)
rvm_all_peaks_std_cpf <- rvm_all_peaks_std[12:22,]
#View(rvm_all_peaks_std_cpf)
rvm_all_peaks_std_d <- rvm_all_peaks_std[23:33,]
#View(rvm_all_peaks_std_d)
rvm_all_peaks_std_treatment <- rvm_all_peaks_std[12:33,]
#View(rvm_all_peaks_std_treatment)



###### not using this at the moment
# convert to long format for ggplot
rvm_all_peaks_gg_wide <- cbind(as.factor(rvm_all_peaks$treatment), rvm_all_peaks_std)
rvm_all_peaks_gg_wide <- as.data.frame(rvm_all_peaks_gg_wide)
#View(rvm_all_peaks_gg_wide)
colnames(rvm_all_peaks_gg_wide)[1] <- "treatment"
rvm_all_peaks_gg_wide$treatment <- as.factor(as.character(rvm_all_peaks_gg_wide$treatment)) 
#View(rvm_all_peaks_gg_wide)
rvm_all_peaks_gg <- melt(rvm_all_peaks_gg_wide, id.vars=c("treatment"))
#View(rvm_all_peaks_gg)

head(rvm_all_peaks_gg)
ggpairs(rvm_all_peaks_gg,                 # Data frame
        aes(color = treatment,  # Color by group (cat. variable)
            alpha = 0.5))     # Transparency

#  ggplot for histogram and metabolite-specific boxplots
head(rvm_all_peaks_gg_wide)
dim(rvm_all_peaks_gg_wide)
ggpairs(
  rvm_all_peaks_gg_wide[,2:28],
  upper = list(continuous = ggally_density, combo = ggally_box_no_facet),
  lower = list(continuous = ggally_points, combo = ggally_dot_no_facet)
)
#######


# the required sample size for covariance matrices depends on the distribution of the variables being measured. 
# For example, if the variables are highly correlated with each other, a smaller sample size may be 
# sufficient to accurately estimate the covariance matrix. However, if the variables are not very 
# correlated, a larger sample size may be needed to accurately capture the relationships among the variables.

#ggm(rvm_all_peaks[,3:14], methods = c("glasso"), community = TRUE,
#    betweenness = TRUE, plot = FALSE, levels = NULL)

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04313-2
# A GGM consists of a network structure of nodes (representing genes, transcripts, metabolites, 
# or proteins), which are inter-connected by edges reflecting significant partial correlations. 
# Partial correlations measure linear associations between pairs of random variables, where the 
# contribution from the remaining variables is adjusted for. Unlike RNs, which are based on 
# Pearson’s correlation, GGMs remove the spurious correlations caused by confounded variables 
# (e.g. when genes share a common regulator). Compared to BNs, GGMs scale up more efficiently 
# to large network analyses and often yield comparable network reconstruction accuracies.

# Partial correlations can be computed from the (standardized) inverse of the covariance matrix 
# (i.e. the precision matrix). In principle, the covariance matrix is unknown and has to be estimated 
# from data. The estimated covariance matrix must be well-conditioned to ensure that its inverse exists, 
# and that numerical (or estimation) errors are not magnified during its computation. 
# The sample covariance, as obtained from a dataset of n samples and p variables, is (i) 
# invertible and well-conditioned when n is greater than p, (ii) invertible but ill-conditioned 
# when n is comparable to p, and (iii) not invertible when n is smaller than p [14]. The last 
# case is known as a ‘high-dimensional problem’, ‘small n, large p’, or ‘n << p This scenario is 
# common in omics’ studies, where often a large set of genes, proteins or metabolites is quantified 
# in relatively few samples.

#all_peaks
# overlapping kernel density plots for each metabolite
# install.packages("GGally")


# paired scatterplots
pairs(rvm_all_peaks_std, col=rvm_all_peaks$treatment)
pairs(rvm_all_peaks_std_control)
pairs(rvm_all_peaks_std_cpf)
pairs(rvm_all_peaks_std_d)

#for everything, we don't have sample sizes to do treatments
n_samples_all <- nrow(rvm_all_peaks_std)
n_samples_control <- nrow(rvm_all_peaks_std_control)
n_samples_cpf <- nrow(rvm_all_peaks_std_cpf)
n_samples_d <- nrow(rvm_all_peaks_std_d)


##### compute correlations covariance
### all together
# very important, correlation is within individuals, not across treatments, so does not
# indicate anything regrarding effects of treatment
all_peaks_cormatrix <- cor_auto(rvm_all_peaks_std)
all_peaks_cov <- cov(rvm_all_peaks_std)
#all_peaks_cormatrix <- var(rvm_all_peaks_std)
#all_peaks_cormatrix <- as.matrix(forceSymmetric(all_peaks_cormatrix))
# check if the matrix is positive definite or not
is.positive.definite(all_peaks_cov, tol=1e-8)
eigen(all_peaks_cov)$values # sll should be positive for positive definite
corrplot(all_peaks_cov, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, diag=F)
# uses glasso
all_peaks_glasso <- qgraph(all_peaks_cov, graph = "glasso", sampleSize = n_samples_all,
                         layout = "spring", refit = FALSE, title = "all_peaks", threshold=T)
plot(all_peaks_glasso)
# Estimate network using glasso needs the covariance matrix
# the glasso function is an implementation of the graphical lasso algorithm, 
# which is a method for estimating sparse precision matrices and the corresponding 
# Gaussian graphical models. The qgraph function is a tool for visualizing 
# and manipulating graphs and networks.
glasso_fit <- glasso(all_peaks_cov, rho=0)
# Extract estimated adjacency matrix from the glasso_fit
adj_matrix <- glasso_fit$wi
# Plot estimated network
all_peaks_qgraph_adjmatrix_cov <- qgraph(adj_matrix, layout = "spring", labels = colnames(all_peaks_cormatrix))
plot(all_peaks_qgraph_adjmatrix_cov)
#library(patchwork)
#(plot(all_peaks_glasso) | plot(all_peaks_glasso))/(plot(all_peaks_qgraph_adjmatrix_cov))
#library(cowplot)
#plot_grid(plot(all_peaks_glasso), plot(all_peaks_glasso), plot(all_peaks_glasso), ncol = 3)

### controls only
# control corrplot works but qgraph FAILS
# control correlation
all_peaks_cormatrix_control <- cor_auto(rvm_all_peaks_std_control)
is.positive.definite(all_peaks_cormatrix_control, tol=1e-8)
eigen(all_peaks_cormatrix_control)$values # sll should be positive for positive definite
corrplot(all_peaks_cormatrix_control, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
# control covariance
all_peaks_cov_control <- cov(rvm_all_peaks_std_control)
is.positive.definite(all_peaks_cov_control, tol=1e-8)
eigen(all_peaks_cov_control)$values # sll should be positive for positive definite
corrplot(all_peaks_cov_control, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

qgraph(all_peaks_cov_control, graph = "glasso", sampleSize = n_samples_control,
       layout = "spring", refit = FALSE, title = "Original")





#combined treatments
all_peaks_cormatrix_treatment <- cor_auto(rvm_all_peaks_std_treatment, forcePD = T)
all_peaks_cormatrix_treatment <- var(rvm_all_peaks_std_treatment)
all_peaks_cormatrix_treatment <- as.matrix(forceSymmetric(all_peaks_cormatrix_treatment))
is.positive.definite(all_peaks_cormatrix_treatment, tol=1e-8)
corrplot(all_peaks_cormatrix_treatment, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(all_peaks_cormatrix_treatment, graph = "glasso", sampleSize = (n_samples_cpf + n_samples_d),
       layout = "spring", refit = FALSE, title = "all_peaks-All treatments")

#chlorpyrifos corrplot and qgraph FAILS
all_peaks_cormatrix_cpf <- cor_auto(rvm_all_peaks_std_cpf) #, forcePD = T
all_peaks_cormatrix_cpf <- var(rvm_all_peaks_std_cpf) #, forcePD = T
all_peaks_cormatrix_cpf <- as.matrix(forceSymmetric(all_peaks_cormatrix_cpf))
is.positive.definite(all_peaks_cormatrix_cpf, tol=1e-8)
corrplot(all_peaks_cormatrix_cpf, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(all_peaks_cormatrix_cpf, graph = "glasso", sampleSize = (n_samples_cpf),
       layout = "spring", refit = FALSE, title = "Original")

# 24d corrplot and qgraph FAILS
all_peaks_cormatrix_d <- cor_auto(rvm_all_peaks_std_d) #, forcePD = T
all_peaks_cormatrix_d <- as.matrix(forceSymmetric(all_peaks_cormatrix_d))
is.positive.definite(all_peaks_cormatrix_d, tol=1e-8)
corrplot(all_peaks_cormatrix_d, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(all_peaks_cormatrix_d, graph = "glasso", sampleSize = (n_samples_d),
       layout = "spring", refit = FALSE, title = "Original")

# trying with 24d and covariance
all_peaks_cov_d <- cov(rvm_all_peaks_std_d)
all_peaks_cov_d <- as.matrix(forceSymmetric(all_peaks_cov_d))
qgraph(all_peaks_cov_d, graph = "glasso", sampleSize = (n_samples_d),
       layout = "spring", refit = FALSE, title = "Original")

# need to implement parameter regularization to overcome the positive definite matrix problem for the individual treatments
# we will use the graphical lasso method (hopefully)
#Friedman, J., Hastie, T., & Tibshirani, R. (2008). 
#Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9(3), 432-441. Chicago
# glasso package
# requires a tuning hyperparameter
# Mazumer, R., Hastie, T., The graphical lasso: New insights and alternatives. Electron J Stat. ; 6: 2125–2149. doi:10.1214/12-EJS740.

#There are several methods to convert a correlation matrix to a positive definite matrix, some of which include:
#Add a small positive constant to the diagonal of the matrix. This is known as "regularization" and is often done with the "Ridge" method of regularization.
#Use a modified Cholesky decomposition, also known as "near-PD" Cholesky decomposition, which uses a regularization term to ensure that the decomposition is positive definite.
#Use the "Matrix square root" method, which calculates the matrix square root of the correlation matrix and then multiplies it by its transpose to obtain a positive definite matrix
#Use the "Eigenvalue adjustment" method, which is to add a small positive constant to the eigenvalues of the correlation matrix and then reconstruct the matrix.
#Use the "Regularized covariance" method, which is to use a shrinkage estimator to estimate the covariance matrix, and then convert it to a correlation matrix.
#It is worth noting that converting a correlation matrix to a positive definite matrix is not always possible, since correlation matrix has constraint that the diagonal elements should be 1, and the off-diagonal elements should be between -1 and 1. Therefore, if the correlation matrix is already positive definite, no conversion is needed.


#EBICglasso from the qgraph pacakge is used to wrap glasso
#https://www.rdocumentation.org/packages/qgraph/versions/1.9.2/topics/EBICglasso

# This function uses the glasso package (Friedman, Hastie and Tibshirani, 2011) to compute a sparse 
# gaussian graphical model with the graphical lasso (Friedman, Hastie \& Tibshirani, 2008). The tuning 
# parameter is chosen using the Extended Bayesian Information criterium (EBIC).
# The glasso is run for 100 values of the tuning parameter logarithmically spaced between the maximal 
# value of the tuning parameter at which all edges are zero, lamba_max, and lambda_max/100. For each of 
# these graphs the EBIC is computed and the graph with the best EBIC is selected. The partial correlation 
# matrix is computed using wi2net and returned. When threshold = TRUE, elements of the inverse variance-covariance 
# matrix are first thresholded using the theoretical bound (Jankova and van de Geer, 2018).
# Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9(3), 432-441. Chicago
# Jerome Friedman, Trevor Hastie and Rob Tibshirani (2011). glasso: Graphical lasso-estimation of Gaussian graphical models. R package version 1.7. http://CRAN.R-project.org/package=glasso
# Foygel, R., & Drton, M. (2010, November). Extended Bayesian Information Criteria for Gaussian Graphical Models. In NIPS (pp. 604-612). Chicago
# Revelle, W. (2014) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston, Illinois, USA, http://CRAN.R-project.org/package=psych Version = 1.4.4.
# Bates, D., and Maechler, M. (2014). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.1-3. http://CRAN.R-project.org/package=Matrix
# Jankova, J., and van de Geer, S. (2018) Inference for high-dimensional graphical models. In: Handbook of graphical models (editors: Drton, M., Maathuis, M., Lauritzen, S., and Wainwright, M.). CRC Press: Boca Raton, Florida, USA.

###everything together


### begin ridge regularization
# defines a regularization term lambda and then add this term to the diagonal of the correlation matrix 
# to obtain a positive definite matrix pd_matrix. The regularization term lambda is a small positive 
# just add a constant to the diagonal of the matrix
epsilon <- 0.0001
all_peaks_cormatrix
all_peaks_cormatrix_ridge <- all_peaks_cormatrix + diag(rep(epsilon, nrow(all_peaks_cormatrix)))
eigen(all_peaks_cormatrix_ridge)$values #its only the upper triangular
is.positive.definite(all_peaks_cormatrix_ridge, tol=1e-8)
### end ridge regularization

### begin cholensky decomposition
# This function checks if the input matrix is positive definite using the chol function. If 
# it is not positive definite, the function will add a small positive constant to the diagonal 
# elements of the matrix until it is positive definite, then return the regularized matrix.
# Function to regularize a correlation matrix
all_peaks_cormatrix_chol <- chol(all_peaks_cormatrix)
is.positive.definite(as.matrix(all_peaks_cormatrix_chol), tol=1e-8)
# make symmetric (machine precision)
all_peaks_cormatrix_chol[lower.tri(all_peaks_cormatrix_chol)] = t(all_peaks_cormatrix_chol)[lower.tri(all_peaks_cormatrix_chol)]
all_peaks_cormatrix_chol
dim(all_peaks_cormatrix_chol)
isSymmetric.matrix(all_peaks_cormatrix_chol)
eigen(all_peaks_cormatrix_chol)$values
is.positive.definite(as.matrix(all_peaks_cormatrix_chol), tol=1e-8)
### end cholensky decomposition

# It's important to note that adding a small positive constant to the diagonal of a correlation matrix 
# to make it positive definite is not a general solution and it's not recommended for all cases, and 
# it's better to consider other options like adding a multiple of the identity matrix, adding a small 
# random noise matrix or using other techniques like the eigenvalue decomposition.

### begin 

all_peaks_glasso <- EBICglasso(all_peaks_cormatrix_chol, n_samples_all, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100, 
                             lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
                             countDiagonal = FALSE, refit = FALSE, threshold = FALSE,
                             verbose = TRUE)
corrplot(all_peaks_glasso, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(all_peaks_glasso, graph = "glasso", sampleSize = n_samples_all,
       layout = "spring", refit = FALSE, title = "Original")

### end everything together