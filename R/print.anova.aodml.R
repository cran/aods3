print.anova.aodml <- function(x, ...) {
  mod <- x$models
  dfr <- x$anova.table
  nam <- rownames(dfr)
  cat("Analysis of Deviance Table\n\n")
  sapply(mod, function(x) cat(x, "\n"))
  cat("\n")
  List <- lapply(dfr, function(x) ifelse(is.na(x), "", format(x, digits = 4)))
  dfr <- as.data.frame(t(do.call("rbind", List)))
  rownames(dfr) <- nam
  print(dfr)
  invisible(dfr)
  }
