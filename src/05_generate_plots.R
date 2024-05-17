## Summary plots see docs/FURTHER_INFO.md for more information.
library(ggplot2)
library(dplyr)
library(ggpubr)
library(GGally)
source("src/util/helper_trt_eff.R")
source("src/util/trt_effect_plots.R")
library(vivid)

## load data without missing covariates imputed
adat <- readRDS(file = "data/analysis_data2.rds")
importance_scores <- readRDS(file = "data/importance.rds")
phi <- readRDS(file = "data/phi.rds")
p_value <- readRDS(file = "data/p_value.rds")

# create a template, users can change it according to their preferences
plot_theme <- theme(
  # Hide panel borders and remove grid lines
  axis.text=element_text(size=12),
  axis.title=element_text(size=14),
  legend.title = element_text(size = 14),
  legend.text = element_text( size = 14),
  strip.text.x = element_text( size = 14),
  plot.title = element_text(size=14))

#####################
## Evidence of TEH ##
#####################
# No plot for visualizing the vidence, just the p-value

######################
## Effect modifiers ##
######################
## variable importance plot
## standard permutation variable importance
p_importance <- ggplot(head(importance_scores %>% arrange(-vi) %>% filter(vi>0), 15),
                       aes(reorder(variable, vi), vi)) +
  geom_col(color="black", fill = "lightblue" ) + 
  xlab("Variables") + ylab("") +
  # ylab("Variable importance\n(average increase in OOB-MSE\n under permutation)")+
  theme_classic() + coord_flip() + plot_theme 

ggsave("reports/figures/var_importance.pdf", p_importance, width = 6, height = 4, units = "in")
# saveRDS(p_importance, file = "reports/figures/var_importance.rds")

## select top 10 variables of interest
vars_to_plot <- rownames(head(importance_scores %>% arrange(-vi), 10))

## scatterplot matrix of variables of interest to identify potential
## correlations/dependencies
p_correlation <- ggpairs(select(adat, vars_to_plot))
ggsave("reports/figures/correlation_graph.pdf", p_correlation, width = 14, height = 14, units = "in")

################################
## Exploratory display of TEH ##
################################
## univariate plot
## - average outcome by treatment group for variables of interest
## - treatment difference for variables of interest
## for paper, choose one continuous and one categorical
vars_to_plot <- c("X1", "X27")
p_uni <- vector(mode = "list", length = 4)
for (i in 1:2) {
  p_uni[[i*2-1]] <- p_yx(
    y.name = "Y", x.name = vars_to_plot[i], trt.name = "trt",
    data = adat, df = 3, family = "gaussian"
  ) + plot_theme
  p_uni[[i*2]] <- p_diffx(
    y.name = "Y", x.name = vars_to_plot[i], trt.name = "trt",
    data = adat, df = 3, family = "gaussian", phi_predictions = NULL
  ) + plot_theme
}

for(i in 1:4){
  ggsave(sprintf("reports/figures/display_uni_%d.pdf", i), p_uni[[i]], width = 6, height = 4, units = "in")
}

########################################################################
## assess interactions using estimated phi's
## note that estimated phis can be outside [-1,1], so that also the
## smoothed phis may be outside [-1,1]

## plot bivariate plots for variables of interest
pairs_to_plot <- list(c("X1", "X9"), c("X1", "X4"))

jj <- 1
for(ii in 1:length(pairs_to_plot)){
  vars_to_plot <- pairs_to_plot[[ii]]
  
  p_biv_1  <- p_yxx(
    y.name = "Y", x.names = vars_to_plot, trt.name = "trt",
    data = adat, df = 2, family = "gaussian"
  ) + plot_theme
  
  ggsave(sprintf("reports/figures/display_bi_%d.pdf", jj), p_biv_1, width = 6, height = 4, units = "in")
  jj <- jj + 1
  
  p_biv_2 <- p_diffxx(
    data = adat, x.names = vars_to_plot, y.name = "Y",
    trt.name = "trt", family = "gaussian", df = 2,
    phi_predictions = NULL, span = 0.8
  ) + plot_theme
  
  ggsave(sprintf("reports/figures/display_bi_%d.pdf", jj), p_biv_2, width = 6, height = 4, units = "in")
  jj <- jj + 1
  
}


## Optional: variable importance & interaction effect plots
## based on the interaction effect calculated from 04_explore_TEH.R
imp_inter <- readRDS("data/importance_inter.rds")
imp_inter_dat <- imp_inter %>%
  tidyr::separate(Variables, sep = "\\*", into = c("var1", "var2")) %>%
  rename(imp = Interaction)
imp_inter_dat2 <- lapply(1:nrow(imp_inter_dat), function(ii){
  xx <- imp_inter_dat[ii, ]
  row2 <- xx
  names(row2) <- c("var2", "var1", "imp")
  rbind(xx, row2)
}) %>% bind_rows()

vars_sel <- unique(c(imp_inter_dat$var1, imp_inter_dat$var2))
imp_dat <- importance_scores %>% filter(variable %in% vars_sel) %>%
  select(var1 = variable, imp = vi) %>% mutate(var2 = var1) %>% rbind(imp_inter_dat2) %>%
  tidyr::pivot_wider(names_from = var2, values_from = imp) %>%
  tibble::column_to_rownames(var = "var1") %>%
  as.matrix()

order_desc <- order(diag(imp_dat), decreasing = T)
imp_order <- imp_dat[order_desc, order_desc]

pdf("reports/figures/var_importance_inter.pdf")
viviHeatmap(mat = imp_order, angle = 45)
dev.off()


