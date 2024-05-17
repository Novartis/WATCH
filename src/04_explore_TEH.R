## Derive variable importance scores and perform global heterogeneity
library(nnls)
library(SuperLearner)
library(permimp)
library(dplyr)
library(clustermq)
library(ggplot2)

adat <- readRDS(file = "data/analysis_regress_imputed.rds")
y <- adat$Y
trt <- select(adat, trt)$trt
trt <- as.numeric(trt == "1")

## First models are build to model the outcome (and for observational
## studies, or when pooling across multiple studies for the treatment
## assignement), based on the stacking procedure as implemented in the
## ['SuperLearner'](https://CRAN.R-project.org/package=SuperLearner) R
## package. Different base models can be utilized for stacking. These
## models are then used to provide pseudo-observations $\phi$ for the
## treatment effect for each patient as in the double-robust learning
## algorithm [Kennedy (2022)](https://arxiv.org/abs/2004.14497). To be
## able to perform statistical inference later the $\phi$ values for
## each patient are created using cross-fitting based on 10-fold
## cross-validation. The prediction for the pseudo-observation for
## every patient is only based on the models where this patient was
## not used for fitting the models (out-of-fold predictions). Note
## that for small data-sets (and or low number of covariates) the
## utilized base models for the super-learner can/should be changed
## (e.g. standard unpenalized regression can also be used).

##############################################################################
## first create pseudo-observations phi (see https://arxiv.org/abs/2004.14497)
##############################################################################

## select covariates to use for modelling the outcome
## make sure to remove outcome
covs_out <- setdiff(colnames(adat), c("Y"))
X_out <- select(adat, all_of(covs_out))
## select covariates to use for modelling treatment (make sure to remove outcome and treatment)
covs_pi <- setdiff(colnames(adat), c("Y", "trt"))
X_pi <- select(adat, all_of(covs_pi))

## super-learners to include
SL.library <- list(
  "SL.glmnet", "SL.xgboost", "SL.ranger"
)


k <- 10
## create k splits, stratified by treatment group
ids <- caret::createFolds(y = adat$trt, k = k)
N <- nrow(adat)
fit_one_fold <- function(i) {
  library(nnls)
  library(SuperLearner)
  train_ids <- setdiff(1:N, ids[[i]])
  test_ids <- ids[[i]]
  ## SuperLearner for outcome
  SL_out <- SuperLearner::SuperLearner(
    Y = y[train_ids],
    X = X_out[train_ids, ],
    SL.library = SL.library, family = "gaussian",
    verbose = FALSE, method = c("method.CC_LS")
  )
  ## SuperLearner for treatment
  SL_pi <- SuperLearner::SuperLearner(
    Y = trt[train_ids],
    X = X_pi[train_ids, ],
    SL.library = SL.library, family = "binomial",
    verbose = FALSE, method = c("method.CC_LS")
  )
  
  ## predictions needed for outcome
  pred_X1 <- pred_X0 <- pred_trt <- X_out[test_ids, ]
  pred_X1$trt <- factor("1", levels = levels(adat$trt))
  pred_X0$trt <- factor("0", levels = levels(adat$trt))
  ## prediction under control
  m0 <- as.numeric(predict(SL_out, newdata = pred_X0)$pred)
  ## prediction under treatment
  m1 <- as.numeric(predict(SL_out, newdata = pred_X1)$pred)
  ## prediction under observed arm
  mtrt <- as.numeric(predict(SL_out, newdata = pred_trt)$pred)
  ## predictions needed for treatment
  pi <- as.numeric(predict(SL_pi, newdata = X_pi[test_ids, ])$pred)
  list(
    i = i, train_ids = train_ids, test_ids = test_ids,
    m0 = m0, m1 = m1, mtrt = mtrt, pi = pi
  )
}
## perform cross-fitting on grid
export <- list(
  adat = adat, X_out = X_out, X_pi = X_pi,
  ids = ids, N = N, SL.library = SL.library, trt = trt, y = y
)
res <- clustermq::Q(fit_one_fold, i = 1:k, n_jobs = k, export = export)

## Calculate pseudo-observations phi
m1 <- m0 <- mtrt <- pi <- numeric(N)
for (i in 1:k) {
  ids <- res[[i]]$test_ids
  m0[ids] <- res[[i]]$m0
  m1[ids] <- res[[i]]$m1
  mtrt[ids] <- res[[i]]$mtrt
  pi[ids] <- res[[i]]$pi
}

## double robust pseudo observations for the treatment difference (see https://arxiv.org/abs/2004.14497)
phi <- m1 - m0 + (trt - pi) / (pi * (1 - pi)) * (y - mtrt)
saveRDS(phi, file = "data/phi.rds")

## Based on the obtained $\phi$ values a global heterogeneity test
## is performed and variable importance is assessed. This is done by
## fitting a random forest based on conditional inference trees (as
## implemented in the ['party'](https://CRAN.R-project.org/package=party)
## R package). One advantage of this approach versus the approach
## implemented in the randomForest package is that conditional inference
## trees do not suffer from variable selection bias (towards variables
## with many split possibilities) when choosing the variables to split
## [Hothorn et al (2006)](https://doi.org/10.1198/106186006X133933).

## Based on the fitted forest then the standard deviation of the
## model-based predictions of the treatment effect for all patients
## are extracted (a large standard deviation would indicate heterogeneity).
## Then the $\phi$ values are permuted against the considered covariates
## and for each permutation the model above is re-fitted and the standard
## deviation of the model-based predictions of the treatment effect is
## calculated. Under the permutation distribution a low standard
## deviation is expected, as the covariates will not explain the $\phi$
## values. A p-value can be extracted based on the proportion of
## permutations standard deviations that are larger than the observed
## standard deviation for the unpermuted data.

## Variable importance is calculated based on the fitted conditional
## random forest using the approaches outlined in [Debeer and Strobl
## (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03622-2),
## as implemented in the
## ['party'](https://CRAN.R-project.org/package=party) and
## ['permimp'](https://CRAN.R-project.org/package=permimp) R
## packages. The basic idea of standard variable importance is to permute
## a covariate and then to observe how much the out-of-bag mean squared
## error for the treatment effect increases (averaged over the trees in
## the forest). To better quantify the unique contribution of a variable
## also conditional permutation variable importance is available, where
## a variable is permuted within strata of variables that are correlated
## with the variable of interest.

########################################################################
## perform global heterogeneity test and obtain variable importance
########################################################################
## Fit the conditional forest model to pseudo observations
n_cov <- ncol(X_pi)
## control parameters for conditional forest
## can assess sensitivity of results to mtry (mtry=5 is the default in cforest,
## n_cov/3 is the default in the randomForest package)
## larger ntree will reduce variability in results (e.g. VI rankings)
control_cforest <- party::cforest_unbiased(mtry = 5, ntree = 500)
fit <- party::cforest(y ~ ., data.frame(X_pi, y = phi), control = control_cforest)
## standard deviation of observed "individual" treatment effects
sd_obs <- sd(predict(fit, OOB = TRUE))
## assess variable importance (to get more stable VI rankings increase nperm)
cf_vi <- permimp::permimp(fit, conditional = FALSE, nperm = 10)
## Optional: The conditional  permutation importance is prohibitive in terms of computational cost.
# cf_vi_cond <- permimp::permimp(fit, conditional = TRUE, nperm = 10)


## assess variability across trees
plot(cf_vi)
plot(cf_vi, type = "box")


## Global heterogeneity test using coin
test <- coin::independence_test(y ~ ., data.frame(X_pi, y = phi), teststat="quadratic")
p_value <- coin::pvalue(test)
p_value

saveRDS(p_value,
        file = "data/p_value.rds"
)

## Optional: can produce p-values for the importance ranking
## p-values for variable importance
## export <- list(X_pi=X_pi, phi=phi, control_cforest = control_cforest, conditional_permutation = TRUE)
## fit_one_perm_cforest <- function(i){
##   library(permimp, lib.loc="/home/bornkbj3/pub_rlib/")
##   phi_perm <- sample(phi)
##   fitperm <- party::cforest(y ~ ., data.frame(X_pi, y=phi_perm), control = control_cforest)
##   cf_vi_perm <- permimp::permimp(fitperm, conditional = conditional_permutation, n_perm=1)
##   cf_vi_perm$values
## }
## res <- clustermq::Q(fit_one_perm_cforest, i=1:n_perm, n_jobs=min(200,n_perm), export=export)
## cf_vi_perm <- do.call("cbind", res)
## p_values <-  numeric(n_cov)
## for(i in 1:n_cov){
##   sm <- sum(cf_vi_perm[i,] > cf_vi$values[i])
##   ## for conditional permutation use cf_vi_cond$values in line above
##   p_values[i] <- (sm+1)/(n_perm+1)
## }
## ## present p-value on "surprise" scale
## data.frame(variable = names(cf_vi$values), surprise = -log2(p_values))

## Optional: can assess univariate association by LR test based on linear model
## pval <- numeric(n_cov)
## lm_fit_null <- lm(y ~ ., data = data.frame(y=phi))
## LR_null <- logLik(lm_fit_null)
## for (j in 1:n_cov) {
##   lm_fit <- lm(y ~ ., data = data.frame(x=X_pi[,j], y=phi))
##   pval[j] <- anova(lm_fit, lm_fit_null)[["Pr(>F)"]][2]
## }
## ## present p-value on "surprise" scale
## data.frame(variable = colnames(X_pi), surprise = -log2(pval))

importance_scores <- data.frame(
  variable = names(cf_vi$values),
  vi = cf_vi$values #, cond_vi = cf_vi_cond$values
)
saveRDS(importance_scores, file = "data/importance.rds")



# 
# ### Optional: interaction variable importance based on partial dependence function
# ### refer paper: https://arxiv.org/pdf/1805.04755.pdf
# ## only based on top nx variables from importance score
nx <- 10
topn <- importance_scores %>% arrange(desc(vi)) %>% slice(1:nx) %>% pull(variable)
comb_var <- combn(topn, 2)

fit_one_inter <- function(i) {
  library(moreparty)
  library(pdp)
  library(vivid)
  GetInteractionStrength(fit, xnames = comb_var[, i])
}

res <- clustermq::Q(fit_one_inter, i = 1:ncol(comb_var), n_jobs = min(200, ncol(comb_var)),
                    export = list(fit = fit, comb_var = comb_var))

saveRDS(res %>% bind_rows(), file = "data/importance_inter.rds")
