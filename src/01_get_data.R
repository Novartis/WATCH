#### generate data from benchtm

library(benchtm)
library(dplyr)

set.seed(202324)

#### data part: data generated on different scenarios saved in "scen_param" from library(benchtm)
data(scen_param)


cases <- scen_param %>% filter(type == "continuous" &
                                 pred == "(X14 > 0.25) & (X1 == 'N')" &
                                 b1_rel == 2)

## generate data using benchtm
dat <- generate_scen_data(scen = cases,  include_truth = F) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at('trt', as.factor)

# set.seed(10)
# cases <- (scen_param %>% filter(type == "continuous"))
# cases[15, ]$pred
# dat <- generate_scen_data(scen = cases[15, ]) %>% 
#   mutate_if(is.character, as.factor) %>%
#   select(-c(trt_effect, prob_diff)) %>%
#   mutate_at('trt', as.factor)


## extract informtion on covariates, response, treatment
X <- dat %>% dplyr::select(starts_with("X"))
Y <- dat$Y
trt <- dat$trt

## randomly generate missing value for baseline covariate with missing proportion 5%, 10%, 20% for X2, X5, X6 
## (for demonstration of data preprocessing)
X2_miss_index <- rbinom(n = nrow(X), size = 1, p = 0.05)
X5_miss_index <- rbinom(n = nrow(X), size = 1, p = 0.10)
X6_miss_index <- rbinom(n = nrow(X), size = 1, p = 0.20)

X[which(X2_miss_index == 1), 2] <- NA
X[which(X5_miss_index == 1), 5] <- NA
X[which(X6_miss_index == 1), 6] <- NA

dat <- cbind(Y, trt, X) %>% mutate_at("trt", as.factor)

# Check is the data folder exists, and if not create it
dir.create(file.path("data"))


saveRDS(dat, file = "data/analysis_data.rds")


