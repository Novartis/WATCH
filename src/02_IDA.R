library(tidyverse)
source('src/util/baseline_plots.R') # util functions for plotting
source('src/util/baseline_dependency.R') # util functions for plotting

## assumes the code in 01_read_data.R has already be used to create data/analysis_data.rds
adat <- readRDS(file = "data/analysis_data.rds")

#######################################################################
## Reproduce existing analyses
#######################################################################
## First step is to cross-check the number of patients in the analysis
## and across treatment groups and reproduce already published analyses
## on the same data-set (e.g. the primary analysis for a specific trial)
## to make sure that the right data are being used in the right way.

## code here

#######################################################################
## Univariate summary of baseline variables
#######################################################################
## Investigate basic properties of the distribution of baseline variables
## (mean, variability/information, skewness, outliers, missingness).
## Results may suggest to use transformations (e.g. log transform) for
## some variables.

# ## first text based summaries
# tab <- tableone::CreateTableOne(data = adat)
# summary(tab)

## histogram/barplots of all variables (save plots as pdf under reports/ for easier review)
## for paper, only draw plot with 2 variables
vars <- c("X5", "X6")
pp_plots <- lapply(vars, plot_density_bar, data = adat)
res <- ggpubr::ggarrange(pp_plots[[1]], pp_plots[[2]], nrow = 1, ncol = 2)

# Check is the reports folder exists, and if not create it
dir.create(file.path("reports"))
# Check is the reports/figures folder exists, and if not create it
dir.create(file.path("reports/figures"))
ggsave("reports/figures/ida_a.pdf",
       res,
       width = 6, height = 4, units = "in"
)

#######################################################################
## Stratified univariate summary of baseline variables
#######################################################################
## Observe baseline summaries stratified by another categorical factor.
## Typical stratification factors are "study or "treatment group".
## In this case we want to compare placebo vs Cosentyx 300mg
## but Cosentyx 300mg is not available in all studies. So both
## treatment and study are of interest as stratification factors

## text based summary
tab <- tableone::CreateTableOne(strata = "trt", data = adat)
print(tab)
# tab <- tableone::CreateTableOne(strata = "STUDYID", data = adat)
# print(tab)

## side-by-side boxplot or stacked bar-plots
vars <- c("X5", "X6")
pp_plots <- lapply(vars, plot_boxbar_by, y.var = "trt", data = adat)

res <- ggpubr::ggarrange(pp_plots[[1]], pp_plots[[2]], nrow = 1, ncol = 2)

ggsave("reports/figures/ida_b.pdf",
       res,
       width = 6, height = 4, units = "in"
)

#######################################################################
## Evalulate missing values/non-informative baseline variables
#######################################################################
## Here missingness & missingness patterns are explored further, in
## addition to assessment of variables with low information (e.g. all
## observations of one variable equal to one value). Both may suggest
## removal of certain variables. In addition for categorical variables
## it is observed whether there are variables with sparsely populated
## categories (may merge those categories with other categories).

## missing variables
p1 <- naniar::gg_miss_var(adat %>% select(X1, X2, X3, X4, X5, X6, X7, X8), 
                          show_pct = TRUE)
## missing variable patterns
p2 <- naniar::vis_miss(adat %>% select(X1, X2, X3, X4, X5, X6, X7, X8)) +
  coord_flip()

res <- ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2)

ggsave("reports/figures/ida_c.pdf",
       res,
       width = 6, height = 4, units = "in"
)



## identification of uninformative variables
nzv <- caret::nearZeroVar(adat[, -1], saveMetrics = TRUE)
head(nzv %>% arrange(desc(nzv)), n = 20)

## identification of categorical variables with sparse categories
## (code in src/util/baseline_dependency.R)
low_freq_categories(adat)

#######################################################################
## Assess dependencies across baseline variables
#######################################################################
## Provides a better understanding of the joint distribution of baseline
## variables and helps identify duplicate (or close-to-duplicate) variables in
## the data. In addition it may help with interpretation of final results.
## To assess dependency between two variables X and Y we calculate
## sqrt(X2/(X2+N)), where X2 is the chi-squared statistic comparing a model for
## p(X|Y) versus p(X).
## If both X and Y are continuous linear models can be used (adjusting for Y linearly)
## and the result is very close to the Pearson correlation.
## If both X and Y are categorical multinomial regression models can be used
## (adjusting for Y as categorical variable). This will give results close to
## Pearson's contingency coefficient.
## For mixed data both multinomial and continuous regression can be used (the
## results assessing p(Y|X) vs p(Y) and p(X|Y) vs p(X) are usually very similar.
## Downside of this approach: For non-continuous data, the maximum achievable
## value can be <1, so assess the results by data-type.

## remove outcome and some problematic variables (e.g. just 1 value)
adat_dep <- adat[, 3:dim(adat)[2]]
dependencies <- get_dep(adat_dep)

## assess dependencies by type of variable comparison
dependencies$results %>%
  filter(Comparison == "continuous" & Correlation > 0.7) %>%
  arrange(desc(Correlation))
## for categorical use lower correlation threshold (as max achievable value can be <1)
dependencies$results %>%
  filter(Comparison == "categorical" & Correlation > 0.5) %>%
  arrange(desc(Correlation))
## mixed comparisons
dependencies$results %>%
  filter(Comparison == "mixed" & Correlation > 0.5) %>%
  arrange(desc(Correlation))

## use hierarchical clustering to assess similarity of variables
hc <- hclust(as.dist(1 - dependencies$cor_mat), method = "average")
pdf(file = "reports/figures/ida_d.pdf", height = 4, width = 6)
plot(hc, hang = -1, xlab = NA, sub = NA, cex = 0.67)
dev.off()

