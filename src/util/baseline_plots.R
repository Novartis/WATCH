## plot spike histogram (continuous vars) or barplot (categorical vars)
plot_density_bar <- function(x.var, data) {
  n_miss <- sum(is.na(data[, x.var]))
  miss_perc <- round(n_miss / dim(data)[1] * 100)
  title <- sprintf("%s, missing: %s (%s%%)", x.var, n_miss, miss_perc)
  data_x <- unlist(data %>% dplyr::select(x = all_of(x.var)) %>% na.omit())
  if (is.numeric(data_x) & length(unique(data_x)) > 5) {
    pp <- data %>%
      dplyr::select(x.var = all_of(x.var)) %>%
      ggplot(mapping = aes(x = (x.var))) +
      geom_histogram(bins = 200) +
      geom_rug()
  } else {
    data_x <- data %>%
      dplyr::select(x = all_of(x.var)) %>%
      na.omit()
    n_miss <- sum(is.na(data[, x.var]))
    range_max <- max(table(data_x))
    pp <- ggplot(data_x, aes(x)) +
      geom_bar() +
      geom_text(stat = "count", aes(label = sprintf("%s (%s %%)", ..count.., round(..count.. / sum(..count..) * 100))), vjust = -0.1) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      ylim(c(0, range_max * 1.1))
  }
  pp +
    theme_bw() +
    xlab(x.var) +
    ylab("") +
    ggtitle(title)
}

### plot boxplot (continuous vars) or barplot (categorical vars) by stratification
plot_boxbar_by <- function(x.var, y.var, data) {
  data_x_y <- data %>%
    dplyr::select(x = all_of(x.var), y = all_of(y.var)) %>%
    na.omit()
  if (unlist(is.numeric(data_x_y$x)) & length(unlist(unique(data_x_y$x))) > 5) { ## continuous
    ## side by side boxplot for x.var using y.var as category
    nns <- table(data_x_y$y)
    pp <- ggplot(data_x_y, aes(x = y, y = x)) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot() +
      geom_jitter(shape = 16, position = position_jitter(0.2), color = "red", alpha = 0.2) +
      xlab("") +
      ylab(x.var) +
      ggtitle(x.var) +
      theme(strip.text.x = element_text(size = 3))
  } else { ## categorical
    pp <- data_x_y %>%
      count(y, x) %>%
      group_by(y) %>%
      mutate(pct = round(prop.table(n) * 100)) %>%
      mutate(n_p = sprintf("%s (%s%%)", n, pct)) %>%
      ggplot(aes(factor(y), n, fill = factor(x))) +
      geom_bar(position = "stack", stat = "identity", alpha = 0.8) +
      geom_text(aes(label = n_p), position = position_stack(vjust = 0.5)) +
      labs(title = x.var, fill = x.var) +
      xlab(x.var) +
      ylab("count") +
      theme(legend.position = "top")
  }
  pp +
    theme_bw()
}

