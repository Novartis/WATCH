## Remove & transform variables (based on IDA results)
## impute remaining missing values in baseline covariates.

library(tidyverse)

adat <- readRDS(file = "data/analysis_data.rds")

#######################################################################
## Variable removal or transformations
#######################################################################
## remove variables due to large missingness
adat <- adat %>% select(-c())
## remove variables due to uninformativeness
adat <- adat %>% select(-c())
## remove duplicated or highly correlated variables
adat <- adat %>% select(-c())
## merge sparse categories
adat <- adat 
## variable transformations
adat <- adat %>%
  mutate(
    log_X5 = log(X5 + 0.01),
    # log_X13 = log(X13),
    log_X13 = log(X13 + 0.01),
    log_X25 = log(X25 + 0.01)
  ) %>% 
  select(-c(X5, X13, X25))


## save transformed data
saveRDS(adat, file = "data/analysis_data2.rds")

#######################################################################
## Impute remaining missing covariates
#######################################################################
## For missing data in the *outcome* variable we propose (for consistency)
## to follow the analytical strategy used in the main pre-specified
## clinical trial analyses for this endpoint, following the decided
## intercurrent event strategies and missing data handling approaches
## (see ICH E9 addendum). In situations, where these approaches are
## complex (e.g. various multiple imputation strategies used) subsequent
## analyses may become time-consuming or infeasible. In this case we
## suggest to use simpler analyses that are still in the spirit of
## the main estimand targeted in the clinical trial analyses.
## This may for example mean using single imputation approaches.

## For missing data in *baseline* variables, we provide following
## considerations and suggestions.
## * Consider to drop covariates for missingness > 10%-20%
## * We don't recommend using multiple imputation for baseline variables
##   as it would result in multiple analysis data-sets and thus more
##   time-consuming subsequent analyses.
## * For single imputation we recommend two methods:
##   (i) Imputation of median or mode of the non-missing values from the
##       same baseline variable.
##   or better
##   (ii) Regression imputation of the missing baseline
##        covariates/biomarkers based on all other variables.
## * While there are different ways of doing (ii) above we
##   recommend to perform multiple imputation to get multiple data-sets
##   with complete baseline variables, but then, rather than using all of 
##   those data-sets in subsequent analyses, taking median/mode on the imputed 
##   values for each missing to get a single dataset. The [`mice`](https://CRAN.R-project.org/package=mice)
##   provides a lot of flexibility and ease of use for this purpose.
## * The imputation of baseline covariates should be independent of
##   outcome or any other post-baseline variables.
## * The missing value indicator method creates extra binary covariates for
##   each covariate with missing data, indicating where a missing value was
##   observed. Running subsequent analyses could then be done using these
##   missing indicators. This allows to assess if the missingness is
##   informative [Groenwold et al (2012)](https://doi.org/10.1503/cmaj.110977).

## We encourage user to run analysis with different imputation methods for a stability check.

## Change variable class to either numeric or factor
adat <- adat %>%
  mutate_if(
    sapply(adat, class) %in% c("integer", "numeric"),
    as.numeric
  ) %>%
  mutate_if(
    sapply(adat, class) %in% c("factor", "character"),
    as.factor
  )


## Missing Value variables check again
missing_value_count <- adat %>%
  summarise(., across(everything(), ~ sum(is.na(.)))) %>% # Use this under R >= 4.0.0
  # summarise_all(funs(sum(is.na(.)))) %>%  # Use this line under R < 4
  as.data.frame() %>%
  `rownames<-`("na_count") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(desc(na_count))

head(missing_value_count)

vars_impute <- c(missing_value_count %>%
                   filter(na_count > 0 & rowname != "Y") %>%
                   select(rowname)) # select covariates includes missing value


## optional step: Add indicator variable for missingness
# adat <- adat %>%
#  mutate_at(vars_impute$rowname, list(missing = ~ is.na(.) * 1.))

## Median/Mode imputation based on same variable
df_univariate_impute <- adat

# Median imputation for numeric variables
for (rn in vars_impute$rowname) {
  if (class(df_univariate_impute[[rn]]) == "numeric") {
    df_univariate_impute[rn][is.na(df_univariate_impute[rn])] <- median(df_univariate_impute[[rn]], na.rm = TRUE)
  }
}

# Mode imputation for factor variables
df_univariate_impute <- randomForest::na.roughfix(df_univariate_impute)

## Single imputation based on regression multiple imputation

# Here multiple imputation of the baseline covariates based on all
# other baseline covariates is performed using the mice package. The
# median (continuous variable) or mode (categorical variable) of the 10
# multiply imputed values are finally used for the single imputed data-set.

## remove outcome and trt (as recommended)
exclude_vars <- c("Y", "trt")
df_regress_impute <- adat
## you can add remove.collinear=FALSE if you have collinear variables in dataset to
## prevent imputed dataset contains NA
imputation <- mice::mice(df_regress_impute %>% dplyr::select(-all_of(exclude_vars)),
                         m = 10, maxit = 5, method = "pmm", printFlag = TRUE, seed = 2020
)
# Find mean/mode of 10 imputations
for (v in vars_impute$rowname) {
  if (class(df_regress_impute[[v]]) == "numeric") {
    imp_median <- apply(imputation$imp[[v]], 1, median)
    df_regress_impute[[v]][as.numeric(names(imp_median))] <- imp_median
  }
  if (class(df_regress_impute[[v]]) == "factor") {
    imp_most <- apply(
      imputation$imp[[v]], 1,
      function(x) {
        names(sort(table(x), decreasing = TRUE)[1])
      }
    )
    df_regress_impute[[v]][as.numeric(names(imp_most))] <- imp_most
  }
}


## Save imputed data
saveRDS(df_univariate_impute, file = "data/analysis_univ_imputed.rds")
saveRDS(df_regress_impute, file = "data/analysis_regress_imputed.rds")

