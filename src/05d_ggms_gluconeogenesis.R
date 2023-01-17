dim(rvm_gluconeogenesis)
View(rvm_gluconeogenesis)

#need to standardize all the numerics
standardize = function(x)
{
  return((x-mean(x))/sd(x))
}

rvm_gluconeogenesis_mat <- rvm_gluconeogenesis[,3:ncol(rvm_gluconeogenesis)]
dim(rvm_gluconeogenesis_mat)
rvm_gluconeogenesis_std_t = apply(rvm_gluconeogenesis_mat,1,standardize)

rvm_gluconeogenesis_std <- t(rvm_gluconeogenesis_std_t)
#View(rvm_gluconeogenesis_std)

rvm_gluconeogenesis_std_control <- rvm_gluconeogenesis_std[1:11,]
rvm_gluconeogenesis_std_cpf <- rvm_gluconeogenesis_std[12:22,]
rvm_gluconeogenesis_std_d <- rvm_gluconeogenesis_std[23:33,]
rvm_gluconeogenesis_std_treatment <- rvm_gluconeogenesis_std[12:33,]

# convert to long format for ggplot
rvm_gluconeogenesis_gg_wide <- cbind(as.factor(rvm_gluconeogenesis$treatment), rvm_gluconeogenesis_std)
colnames(rvm_gluconeogenesis_gg_wide)[1] <- "treatment"
rvm_gluconeogenesis_gg <- melt(rvm_gluconeogenesis_gg_wide, id.vars=c("treatment"))
#View(rvm_gluconeogenesis_gg)



# the required sample size for covariance matrices depends on the distribution of the variables being measured. 
# For example, if the variables are highly correlated with each other, a smaller sample size may be 
# sufficient to accurately estimate the covariance matrix. However, if the variables are not very 
# correlated, a larger sample size may be needed to accurately capture the relationships among the variables.

#ggm(rvm_gluconeogenesis[,3:14], methods = c("glasso"), community = TRUE,
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

#gluconeogenesis
# overlapping kernel density plots for each metabolite
# install.packages("GGally")
library(GGally)

ggpairs(rvm_gluconeogenesis_gg,                 # Data frame
        aes(color = treatment,  # Color by group (cat. variable)
            alpha = 0.5))     # Transparency

# paired scatterplots
pairs(rvm_gluconeogenesis_std, col=rvm_gluconeogenesis$treatment)
#for everything, we don't have sample sizes to do treatments
n_samples_all <- nrow(rvm_gluconeogenesis_std)
n_samples_control <- nrow(rvm_gluconeogenesis_std_control)
n_samples_cpf <- nrow(rvm_gluconeogenesis_std_cpf)
n_samples_d <- nrow(rvm_gluconeogenesis_std_d)

library(matrixcalc)
library(Matrix)
# compute correlations covariance
gluconeogenesis_cormatrix <- cor_auto(rvm_gluconeogenesis_std)
gluconeogenesis_cormatrix <- as.matrix(forceSymmetric(gluconeogenesis_cormatrix))
is.positive.definite(gluconeogenesis_cormatrix, tol=1e-8)
# very important, correlation is within individuals, not across treatments, so does not
# indicate anything regrarding effects of treatment
corrplot(gluconeogenesis_cormatrix, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(gluconeogenesis_cormatrix, graph = "glasso", sampleSize = n_samples_all,
       layout = "spring", refit = FALSE, title = "gluconeogenesis")

# control corrplot works but qgraph FAILS
gluconeogenesis_cormatrix_control <- cor_auto(rvm_gluconeogenesis_std_control)
gluconeogenesis_cormatrix_control <- as.matrix(forceSymmetric(gluconeogenesis_cormatrix_control))
is.positive.definite(gluconeogenesis_cormatrix_control, tol=1e-8)
corrplot(gluconeogenesis_cormatrix_control, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(gluconeogenesis_cormatrix_control, graph = "glasso", sampleSize = n_samples_control,
       layout = "spring", refit = FALSE, title = "Original")

#combined treatments
gluconeogenesis_cormatrix_treatment <- cor_auto(rvm_gluconeogenesis_std_treatment, forcePD = T)
gluconeogenesis_cormatrix_treatment <- as.matrix(forceSymmetric(gluconeogenesis_cormatrix_treatment))
is.positive.definite(gluconeogenesis_cormatrix_treatment, tol=1e-8)
corrplot(gluconeogenesis_cormatrix_treatment, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(gluconeogenesis_cormatrix_treatment, graph = "glasso", sampleSize = (n_samples_cpf + n_samples_d),
       layout = "spring", refit = FALSE, title = "Original")

#chlorpyrifos corrplot and qgraph FAILS
gluconeogenesis_cormatrix_cpf <- cor_auto(rvm_gluconeogenesis_std_cpf) #, forcePD = T
gluconeogenesis_cormatrix_cpf <- as.matrix(forceSymmetric(gluconeogenesis_cormatrix_cpf))
is.positive.definite(gluconeogenesis_cormatrix_cpf, tol=1e-8)
corrplot(gluconeogenesis_cormatrix_cpf, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(gluconeogenesis_cormatrix_cpf, graph = "glasso", sampleSize = (n_samples_cpf),
       layout = "spring", refit = FALSE, title = "Original")

# 24d corrplot and qgraph FAILS
gluconeogenesis_cormatrix_d <- cor_auto(rvm_gluconeogenesis_std_d) #, forcePD = T
gluconeogenesis_cormatrix_d <- as.matrix(forceSymmetric(gluconeogenesis_cormatrix_d))
is.positive.definite(gluconeogenesis_cormatrix_d, tol=1e-8)
corrplot(gluconeogenesis_cormatrix_d, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(gluconeogenesis_cormatrix_d, graph = "glasso", sampleSize = (n_samples_d),
       layout = "spring", refit = FALSE, title = "Original")

# trying with 24d and covariance
gluconeogenesis_cov_d <- cov(rvm_gluconeogenesis_std_d)
gluconeogenesis_cov_d <- as.matrix(forceSymmetric(gluconeogenesis_cov_d))
qgraph(gluconeogenesis_cov_d, graph = "glasso", sampleSize = (n_samples_d),
       layout = "spring", refit = FALSE, title = "Original")

# need to implement parameter regularization to overcome the positive definite matrix problem for the individual treatments
# we will use the graphical lasso method (hopefully)
#Friedman, J., Hastie, T., & Tibshirani, R. (2008). 
#Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9(3), 432-441. Chicago
# glasso package
# requires a tuning hyperparameter
# Mazumer, R., Hastie, T., The graphical lasso: New insights and alternatives. Electron J Stat. ; 6: 2125–2149. doi:10.1214/12-EJS740.

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
gluconeogenesis_glasso <- EBICglasso(gluconeogenesis_cormatrix, n_samples_all, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100, 
                          lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
                          countDiagonal = FALSE, refit = FALSE, threshold = FALSE,
                          verbose = TRUE)
corrplot(gluconeogenesis_glasso, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
qgraph(gluconeogenesis_glasso, graph = "glasso", sampleSize = n_samples_all,
       layout = "spring", refit = FALSE, title = "Original")


gluconeogenesis_glasso_control <- EBICglasso(gluconeogenesis_cormatrix_control, n_samples_control, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100, 
                                  lambda.min.ratio = 20, returnAllResults = FALSE, checkPD = TRUE, 
                                  countDiagonal = FALSE, refit = FALSE, threshold = FALSE,
                                  verbose = TRUE)
corrplot(gluconeogenesis_glasso_control, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

gluconeogenesis_glasso_cpf <- EBICglasso(gluconeogenesis_cormatrix_cpf, n_samples_cpf, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100, 
                              lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
                              countDiagonal = FALSE, refit = FALSE, threshold = FALSE,
                              verbose = TRUE)
corrplot(gluconeogenesis_glasso_cpf, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

gluconeogenesis_glasso_d <- EBICglasso(gluconeogenesis_cormatrix_d, n_samples_d, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100, 
                            lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
                            countDiagonal = FALSE, refit = FALSE, threshold = FALSE,
                            verbose = TRUE)
corrplot(gluconeogenesis_glasso_d, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
