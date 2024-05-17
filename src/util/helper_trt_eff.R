## return a plot, with smoothing curve (y versus x) as well as histogram on bottom (based on x)
## if x categorical, return bar plot with count
#' Plot (smoothed) treatment difference against baseline covariate
#'
#' This function plots the (smoothed) treatment difference (in average
#' outcomes) against a baseline covariate. For numeric baseline
#' covariates a regression spline is used (with user-specified df). For
#' categorical baseline covariates just means are presented.
#'
#' @param data output from p_diffx_cal, include columns x, trt_eff, lower, upper 
#' (lower and upper bound for treatment effect), lev(optional), pred(optional): smoothed phi
#' @param dat_raw include x with lev variable
#' @param x.name Name of baseline variable in data
#' @param y.name Name of treatment variable in data (code assumes
#'   there are only two treatments)
#' @param lev.name facet variable name
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @param span Span to use for smoothing the phis
#' @export
#' @examples

comb_plot <- function(data, dat_raw, x.name, y.name, lev.name, 
                      diff_lab = "probability diff"){
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  dat <- data 
  
  if(is.null(lev.name)){
    title_txt <- paste0(diff_lab, " of ", y.name, " on ", x.name)
  }else{
    title_txt <- paste0(diff_lab, " of ", y.name, " on ", x.name, " by ", lev.name)
  }
  
  
  if(is.numeric(dat$x)){
    p_effect <- ggplot(aes(x, trt_eff), data = dat) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "longdash") +
      common_theme +
      labs(
        # title = title_txt,
        # subtitle = "Average with 95% CI",
        x = x.name, y = diff_lab)
    if(!is.null(dat$phi_pred)){
      p_effect <- p_effect + geom_line(aes(x, phi_pred), lty = 2, alpha = 0.5)
    }
    
    p_dist <- ggplot(dat_raw, aes(x = x)) +
      geom_histogram() +
      theme_bw() +
      xlab(x.name)
    
    if(!is.null(lev.name)){
      dat$lev <- data[, "lev"]
      p_effect <- p_effect + facet_grid(~lev)
      p_dist <- p_dist + facet_grid(~x2)
    }
    p_plot <- cowplot::plot_grid(p_effect +
                                   theme(
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank()
                                   ),
                                 p_dist +
                                   theme(strip.background = element_blank(), strip.text.x = element_blank()),
                                 nrow = 2, rel_heights = c(3, 1), align = "v", axis = "lr"
    )
  }else{
    # browser()
    p_plot <- dat %>%
      mutate( x = factor(x, levels = unique(dat$x))) %>%
      ggplot(aes(x, trt_eff)) +
      geom_point() +
      theme_bw() +
      geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5, width = 0.1) +
      geom_hline(yintercept = 0, linetype = "longdash") +
      labs(
        # title = title_txt,
        # subtitle = paste0("Average with 95% CI"),
        x = x.name, y = diff_lab)
    
    if(!is.null(dat$phi_pred)){
      p_plot <- p_plot + geom_point(aes(x, phi_pred), pch = 4, alpha = 0.5, data = dat)
    }
    
    if(!is.null(lev.name)){
      dat$lev <- data[, "lev"]
      p_plot <- p_plot + facet_grid(~lev, scale = "free")
    }
  }
  return(p_plot)
}


## functions to help calculate the dataframe used for p_diffx and p_diffxx
#'
#' This function calculate treatment difference 
#' (response difference or loglogs ratio difference) given x and trt
#'
#' @param dat a data frame with variable x, trt, and y
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @param phi_predictions Predicted phi's obtained from the DR learner (optional)
#' @param span Span to use for smoothing the phis
#' @export
#' @examples
p_diffx_cal <- function(dat,
                        df,
                        family,
                        span) {
  levels_trt <- levels(dat$trt)
  if (is.numeric(dat$x)) {
    ## code based on glm and contrast
    x_seq <-
      seq(min(dat$x, na.rm = T), max(dat$x, na.rm = T), length.out = 100)
    glm_fit1 <-
      glm(
        y ~ splines::ns(x, df = df),
        data = filter(dat, trt == levels_trt[1]),
        family = family
      )
    glm_fit2 <-
      glm(
        y ~ splines::ns(x, df = df),
        data = filter(dat, trt == levels_trt[2]),
        family = family
      )
    pred_dat <- data.frame(x = x_seq)
    pred1 <-
      predict(glm_fit1, pred_dat, type = "response", se.fit = TRUE)
    pred2 <-
      predict(glm_fit2, pred_dat, type = "response", se.fit = TRUE)
    plot_dat <-
      data.frame(
        x = x_seq,
        trt_eff = pred1$fit - pred2$fit,
        se = sqrt(pred1$se ^ 2 + pred2$se ^ 2)
      ) %>%
      mutate(lower = trt_eff - 2 * se, upper = trt_eff + 2 * se)
  } else {
    n_count <- dat %>%
      mutate( x = factor(x, levels = unique(dat$x))) %>%
      group_by(x) %>%
      summarise(n = n()) %>%
      mutate(x_label = paste0(x, "(n = ", n, ")")) %>%
      ungroup()
    
    dat <- dat %>%
      mutate(x = factor(x, levels = n_count$x, labels = n_count$x_label))
    
    x_lev <- n_count$x_label
    
    glm_fit <- glm(y ~ trt * x, data = dat, family = family)
    preds <- pred_se <- numeric(length(x_lev))
    for (i in 1:length(x_lev)) {
      glm_fit <-
        glm(y ~ trt,
            data = filter(dat, x == x_lev[i]),
            family = family)
      pred1 <- predict(
        glm_fit,
        newdata = data.frame(trt = levels_trt[1], x = x_lev[i]),
        se.fit = TRUE,
        type = "response"
      )
      pred2 <- predict(
        glm_fit,
        newdata = data.frame(trt = levels_trt[2], x = x_lev[i]),
        se.fit = TRUE,
        type = "response"
      )
      preds[i] <- pred1$fit - pred2$fit
      pred_se[i] <- sqrt(pred1$se.fit ^ 2 + pred2$se.fit ^ 2)
    }
    plot_dat <-
      data.frame(
        x = factor(x_lev, levels  = x_lev),
        trt_eff = preds,
        se = pred_se
      ) %>%
      mutate(lower = trt_eff - 2 * se, upper = trt_eff + 2 * se)
  }
  
  if (is.numeric(dat$x) & (!is.null(dat$phi))) {
    ## refit phi using smoothing curve
    fit <- loess(phi ~ x, data = dat, span = span)
    plot_dat$phi_pred <- predict(fit, newdata = plot_dat)
  }else if(!is.numeric(dat$x) & (!is.null(dat$phi))){
    mns <- tapply(dat$phi, dat$x, mean)
    pdat2 <- data.frame(x = names(mns), phi_pred = mns)
    plot_dat <- full_join(pdat2, plot_dat, by = "x")
  }
  
  return(plot_dat)
}




comb_plot_yx <- function(data, x.name, y.name, trt.name, lev.name, family, df = 2){
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  dat <- data 
  
  if(is.null(lev.name)){
    title_txt <- paste0(y.name, " on ", x.name)
  }else{
    title_txt <- paste0(y.name, " on ", x.name, " by ", lev.name)
  }
  
  
  if(is.numeric(dat$x)){
    
    p_effect <- dat %>%
      ggplot(aes(x, y, col = trt)) +
      geom_point() +
      geom_smooth(
        method = "glm", method.args = list(family = family),
        formula = y ~ splines::ns(x, df = df)
      ) +
      common_theme +
      labs(
        # title = title_txt,
        # subtitle = "Average with 95% CI",
        x = x.name, y = y.name, color = trt.name
      )
    p_dist <- ggplot(dat, aes(x = x)) +
      geom_histogram() +
      theme_bw() +
      xlab(x.name) 
    
    if(!is.null(lev.name)){
      dat$lev <- dat[, "lev"]
      p_effect <- p_effect + facet_grid(~lev)
      p_dist <- p_dist + facet_grid(~lev)
    }
    p_plot <- cowplot::plot_grid(p_effect +
                                   theme(
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank()
                                   ),
                                 p_dist +
                                   theme(strip.background = element_blank(), strip.text.x = element_blank()),
                                 nrow = 2, rel_heights = c(3, 1), align = "v", axis = "lr"
    )
  }else{
    # browser()
    
    p_plot <- dat %>%
      ggplot(aes(x, y, col = trt, group = trt)) +
      stat_summary(fun = mean, geom = "point") +
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1) +
      common_theme +
      labs(
        # title = title_txt,
        # subtitle = "Average with 95% CI",
        x = x.name, y = y.name, color = trt.name
      )
    
    if(!is.null(lev.name)){
      dat$lev <- data[, "lev"]
      p_plot <- p_plot + facet_grid(~lev, scale = "free")
    }
  }
  return(p_plot)
}

