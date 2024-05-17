get_cor <- function(x, y, type_x, type_y) {
  ## function to calculate the correlation between variables of mixed
  ## type. In case both x,y are continuous, the standard correlation
  ## is returned. Otherwise an alternative that would reduce to
  ## Pearson's correlation for continuous data
  Pearson_C <- function(X2, N) {
    X2 <- max(X2, 0)
    sqrt(X2 / (N + X2))
  }
  N <- length(x)
  type_x_cont <- type_x == "continuous"
  type_y_cont <- type_y == "continuous"
  if (type_x_cont & type_y_cont) {
    fit <- lm(y ~ x)
    fit0 <- lm(y ~ 1)
    X2 <- as.numeric(2 * (logLik(fit) - logLik(fit0)))
    return(Pearson_C(X2, N)) # should roughly equal cor(x,y)
  }
  if (!type_x_cont & !type_y_cont) { # both categorical
    fit <- nnet::multinom(y ~ x, trace = FALSE)
    fit0 <- nnet::multinom(y ~ 1, trace = FALSE)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    return(Pearson_C(X2, N))
  }
  if (type_x_cont & !type_y_cont) {
    fit <- nnet::multinom(y ~ x, trace = FALSE) ## multinomial regression
    fit0 <- nnet::multinom(y ~ 1, trace = FALSE)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P1 <- Pearson_C(X2, N)
    fit <- lm(x ~ y)
    fit0 <- lm(x ~ 1)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P2 <- Pearson_C(X2, N)
    ## cat(sprintf("P1: %.2f, P2: %.2f\n", P1, P2))
    return((P1 + P2) * 0.5)
  }
  if (!type_x_cont & type_y_cont) {
    fit <- nnet::multinom(x ~ y, trace = FALSE) ## multinomial regression
    fit0 <- nnet::multinom(x ~ 1, trace = FALSE)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P1 <- Pearson_C(X2, N)
    fit <- lm(y ~ x)
    fit0 <- lm(y ~ 1)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P2 <- Pearson_C(X2, N)
    ## cat(sprintf("P1: %.2f, P2: %.2f\n", P1, P2))
    return((P1 + P2) * 0.5)
  }
}

get_dep <- function(dat) {
  var_type <- sapply(dat, function(x) {
    cl <- class(x)
    if (cl %in% c("factor", "ordered", "logical", "character")) {
      return("categorical")
    }
    if (cl %in% c("numeric", "integer", "double")) {
      return("continuous")
    }
  })
  n_vars <- ncol(dat)
  out <- matrix(0, nrow = n_vars, ncol = n_vars)
  res <- vector("list", 0.5 * n_vars * (n_vars - 1))
  nams <- names(dat)
  colnames(out) <- rownames(out) <- names(dat)
  z <- 1
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      tmp <- get_cor(dat[, i, drop = TRUE], dat[, j, drop = TRUE], var_type[i], var_type[j])
      if (inherits(tmp, "try-error")) {
        txt <- sprintf("Problem for variables: %s and %s\n", nams[i], nams[j])
        cat(txt, attr(tmp, "condition")$message, "\n")
        cat("Contingency table for variables looks like this (missing values removed)\n")
        print(table(dat[, i], dat[, j]))
        break
      } else {
        out[i, j] <- out[j, i] <- tmp
        if (var_type[i] == var_type[j]) {
          comparison <- var_type[i]
        } else {
          comparison <- "mixed"
        }
        res[[z]] <- c(nams[i], nams[j], tmp, comparison)
        z <- z + 1
      }
    }
  }
  res <- as.data.frame(do.call("rbind", res))
  names(res) <- c("Variable 1", "Variable 2", "Correlation", "Comparison")
  res[, 3] <- as.numeric(res[, 3])
  list(cor_mat = out, results = res)
}

low_freq_categories <- function(dat, perc) {
  var_type <- sapply(dat, function(x) {
    cl <- class(x)
    if (cl %in% c("factor", "ordered", "logical", "character")) {
      return("categorical")
    }
    if (cl %in% c("numeric", "integer", "double")) {
      return("continuous")
    }
  })
  ind <- var_type == "categorical"
  vars <- sapply(dat[, ind], function(x) {
    tab <- table(x)
    freq_tab <- tab / sum(tab)
    min(freq_tab)
  })
  vars <- vars[order(vars, decreasing = FALSE)]
  data.frame(var = names(vars), size_sparsest_cat = as.numeric(vars))
}

