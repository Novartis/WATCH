##########################
## functions to generate plot for response by treatment/treatment effect v.s. covariates(one or two variables)
## to explore the relationship
## Authors: sophie.sun@novartis.com and bjoern.bornkamp@novartis.com
############################

#' Plot (smoothed) average outcome against baseline covariate by
#' treatment group
#'
#' This function plots the (smoothed) average outcome against baseline
#' covariate by treatment group. For numeric baseline covariates a
#' regression spline is used (with user-specified df). For categorical
#' baseline covariates just means are presented.
#'
#' @param y.name Name of outcome variable in data
#' @param x.name Name of baseline variable in data
#' @param trt.name Name of treatment variable in data
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @export
#' @examples
p_yx <- function(y.name, x.name, trt.name, data, df = 2, family) {
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  # browser()
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.name))) %>% na.omit()
  colnames(dat) <- c("y", "trt", "x")
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(!is.numeric(dat$x)){
    n_count <- dat %>%
      mutate( x = factor(x, levels = unique(dat$x))) %>%
      group_by(x) %>%
      summarise(n = n()) %>%
      mutate(x_label = paste0(x, "(n = ", n, ")")) %>%
      ungroup()
    
    dat <- dat %>%
      mutate(x = factor(x, levels = n_count$x, labels = n_count$x_label))
  }
  
  comb_plot_yx(data = dat, x.name, y.name, trt.name, lev.name = NULL, family = family, df = df)
}


#' Plot (smoothed) treatment difference against baseline covariate
#'
#' This function plots the (smoothed) treatment difference (in average
#' outcomes) against a baseline covariate. For numeric baseline
#' covariates a regression spline is used (with user-specified df). For
#' categorical baseline covariates just means are presented.
#'
#' @param y.name Name of outcome variable in data
#' @param x.name Name of baseline variable in data
#' @param trt.name Name of treatment variable in data (code assumes
#'   there are only two treatments)
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @param facet_title Include a facet (useful for p_trt_effect_xx function)
#' @param label_y Include the label of the y axis (useful for p_trt_effect_xx function)
#' @export
#' @examples
p_diffx <- function(y.name, x.name, trt.name, data, family, df = 3,
                    phi_predictions = NULL, span = 0.75) {
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  # browser()
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.name)))
  colnames(dat) <- c("y", "trt", "x")
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(!is.null(phi_predictions)){
    dat$phi <- phi_predictions
  }
  plot_dat <- p_diffx_cal(dat, df, family, span = span)
  
  
  comb_plot(data = plot_dat, dat_raw = dat %>% select(x),
            x.name = x.name, y.name = y.name, lev.name = NULL, 
            diff_lab = "Probability difference")           
}



###############################################################################
p_diffxx <- function(data, x.names, y.name, trt.name, family, df = 3,
                     phi_predictions = NULL, span = 0.75) {
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  # browser()
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.names))) 
  colnames(dat) <- c("y", "trt", "x1", "x2")
  
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(!is.null(phi_predictions)){
    dat$phi <- phi_predictions
  }
  
  dat <- dat %>% drop_na()
  
  if(is.numeric(dat$x1)){
    ## x1 continuous, x2 continuous, categorize x2 and then facet by x2, facet_name is the original name of x2
    dat = as.data.frame(dat) %>% rename(x = x1)
    if(is.numeric(dat$x2)){
      dat[,'x2'] = cut( dat[,'x2'] , breaks=quantile( dat[,'x2'], seq(0.00,1,1/3), na.rm = TRUE),  include.lowest = TRUE)
    }
    
    n_count <- dat %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[2], ": ", x2, " (n = ", n, ")"))
    
    dat <- dat %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label))
    
    dat_plot <- lapply(split(dat, dat[, "x2"]), function(xx_dat){
      p_diffx_cal(dat = xx_dat, df = df, family = family,
                  span = span) %>% mutate(lev = xx_dat[1, "x2"])
    }) %>% bind_rows()
    
    comb_plot(data = dat_plot, dat_raw = dat %>% select(x, x2),
              x.name = x.names[1], y.name = y.name, lev.name = x.names[2], 
              diff_lab = "Probability difference")
    
  }else if (!is.numeric(dat$x1) & is.numeric(dat$x2)){
    ## x1 categorical, x2 continuous, facet by x1, facet_name is the original name of x1
    dat <- as.data.frame(dat) %>% rename(x = x2, x2 = x1)
    ## x continuous, x2 categorical
    
    n_count <- dat %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[1], ": ", x2, " (n = ", n, ")"))
    dat <- dat %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label))
    
    dat_plot <- lapply(split(dat, dat[, "x2"]), function(xx_dat){
      p_diffx_cal(dat = xx_dat, df = df, family = family,
                  span = span) %>% mutate(lev = xx_dat[1, "x2"])
    }) %>% bind_rows()
    
    comb_plot(data = dat_plot, dat_raw = dat %>% select(x, x2),
              x.name = x.names[2], y.name = y.name, lev.name = x.names[1], 
              diff_lab = "Probability difference")
    
  } else{
    dat <- as.data.frame(dat) %>% rename(x = x1)
    
    n_count <- dat %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[2], ": ", x2, " (n = ", n, ")"))
    dat <- dat %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label))
    
    ## x categorical, x2 categorical
    dat_plot <- lapply(split(dat, dat[, "x2"]), function(xx_dat){
      p_diffx_cal(dat = xx_dat, df = df, family = family,
                  span = span) %>% mutate(lev = xx_dat[1, "x2"])
    }) %>% bind_rows() %>% mutate(x = as.character(x)) %>% arrange(lev, x)
    
    comb_plot(data = dat_plot, dat_raw = dat %>% select(x, x2),
              x.name = x.names[1], y.name = y.name, lev.name = x.names[2], 
              diff_lab = "Probability difference")
  } 
}


## bivariate outcome plots
## combination of two categorical variables or one continuous variable and one categorical variable

#' Plot (smoothed) treatment difference against two baseline 
#' covariate, to assess potential interactions
#'
#' This function plots the estimated outcome estimated from original outcome against two
#' baseline covariates. If both baseline variables are categorical, the phi's are smoothed using loess
#' and displayed using a contour plotit is simple average with CI. If one baseline variable is continuous and one categorical
#' a regression spline is used (with user-specified df) for each category of the second variable.
#' 
#'
#' @param y.name Name of outcome variable in data
#' @param x.names Names of baseline variable in data
#' @param trt.name Name of treatment variable in data (code assumes
#'   there are only two treatments)
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @export
#' @examples
p_yxx <- function(y.name, x.names, trt.name, data, df = 2, family = "binomial") {
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.names))) %>% drop_na()
  colnames(dat) <- c("y", "trt", "x1", "x2")
  
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(is.numeric(dat$x1)){
    ## x1 continuous, x2 continuous, categorize x2 and then facet by x2, facet_name is the original name of x2
    dat_new <- as.data.frame(dat) %>% rename(x = x1)
    if(is.numeric(dat$x2)){
      dat_new[,'x2'] = cut( dat_new[,'x2'] , breaks=quantile( dat_new[,'x2'], seq(0.00,1,1/3), na.rm = TRUE),  include.lowest = TRUE)
    }
    
    n_count <- dat_new %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[2], ": ", x2, " (n = ", n, ")"))
    
    dat_new <- dat_new %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label)) %>% rename(lev = x2)
    
    comb_plot_yx(data = dat_new, x.name = x.names[1], y.name, trt.name, 
                 lev.name = x.names[2], family = family, df = df)
    
  }else if (!is.numeric(dat$x1) & is.numeric(dat$x2)){
    ## x1 categorical, x2 continuous, facet by x1, facet_name is the original name of x1
    dat_new <- as.data.frame(dat) %>% rename(x = x2, x2 = x1)
    ## x continuous, x2 categorical
    
    n_count <- dat_new %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[1], ": ", x2, " (n = ", n, ")"))
    dat_new <- dat_new %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label)) %>% rename(lev = x2)
    
    comb_plot_yx(data = dat_new, x.name = x.names[2], y.name, trt.name, 
                 lev.name = x.names[1], family = family, df = df)
    
  } else{
    dat_new <- as.data.frame(dat) %>% rename(x = x1)
    
    dat_new <- lapply(split(dat_new, dat_new$x2), function(xx){
      n_count <- xx %>%
        group_by(x) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        mutate(x_label = paste0(x, " (n = ", n, ")"))
      xx %>% mutate(x = factor(x, levels = n_count$x, labels = n_count$x_label)) %>% 
        mutate(x2 = paste0(x.names[2], ": ", xx$x2[1], " (n = ", nrow(xx), ")"))
    }) %>% bind_rows() %>% rename(lev = x2)
    
    comb_plot_yx(data = dat_new, x.name = x.names[1], y.name, trt.name, 
                 lev.name = x.names[2], family = family, df = df)
  } 
}


########################################################################
## helper functions to create exhaustive subgroup plot (funnel plot)


#' function to guess variable type (numeric or categorical)
#'
#' @param x Data frame
#' @export
#' @examples
guess_type <- function(x) {
  N <- length(x)
  uN <- length(unique(x))
  ind <- !is.na(x) ## remove missings
  oldw <- getOption("warn")
  options(warn = -1)
  xx <- as.double(x[ind])
  options(warn = oldw)
  mn <- mean(is.na(xx))
  if (((uN > 8) | (N < 40 & uN / N > 0.4)) & mn < 0.01) { ## assume it's numeric
    return("numeric")
  } else {
    if (uN < 40) {
      return("categorical")
    } else {
      return("undigestable")
    }
  }
}

#' function to categorize a continuous variable in X categorizes
#' according to quantiles and returns in a list (i) the categorized
#' variable (ii) the category labels and (iii) character variable with
#' describing the link between labels and cut-offs
#'
#' @param x Data frame
#' @param nr Number of categories to categorize a numerical variable
#' @param var_nam Intended name for the variable
#' @export
#' @examples

categ_var <- function(x, nr, var_nam) {
  psq <- seq(0, 1, length = nr + 1)
  qq <- quantile(x, psq, na.rm = TRUE)
  delta <- diff(range(x, na.rm = TRUE))
  qq[1] <- min(x, na.rm = TRUE) - 0.02 * delta
  qq[nr + 1] <- max(x, na.rm = TRUE) + 0.02 * delta
  brks <- unique(qq)
  if(length(brks) == 3){
    labs <- c("low", "high")
    char <- sprintf("%s: low (<= %s) | high (> %s)",
                    var_nam, brks[2], brks[2])
  }
  if(length(brks) == 4){
    labs <- c("low", "mid", "high")
    char <- sprintf("%s: low (<= %s) | mid (%s, %s] | high (> %s)",
                    var_nam, brks[2], brks[2], brks[3], brks[3])
  }
  ct <- cut(x, brks, labels = labs)
  list(ct, labs, char)
}

#' function to produce categorical subgroup variables
#'
#' @param data Data frame
#' @export
#' @examples
get_subg_var <- function(data) {
  types <- sapply(data, guess_type)
  nams <- names(types)
  ind <- types == "undigestable"
  if (any(ind)) {
    types <- types[-ind]
    nams <- nams[-ind]
  }
  lst <- levs <- vector("list", length = length(nams))
  chars <- character(length(nams))
  for (i in 1:length(nams)) {
    if (types[i] == "numeric") 
      tmp <- categ_var(data[, get(nams[i])], 3, nams[i])
    if (types[i] == "categorical"){
      labs <- as.character(unique(data[, get(nams[i])]))
      tmp <- list(data[, get(nams[i])],
                  labs,
                  sprintf("%s: %s", nams[i], paste0(labs, collapse=" | ")))
    }
    lst[[i]] <- tmp[[1]]
    levs[[i]] <- tmp[[2]]
    chars[i] <- tmp[[3]]
  }
  names(lst) <- nams
  names(levs) <- nams
  list(subgr_vars = lst, subgr_levels = levs, chars = chars)
}
