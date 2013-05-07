anova.aodml <- function(object, ...) {

  #### copied from anova method for lmer (packages Matrix / lme4)

  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- sapply(dots, inherits, "aodml")
  mods <- c(list(object), dots[modp])
  names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)

	####
	
	n <- length(mods)
  if(n < 2)
    stop("At least 2 valid models are needed.")
  
	dfr <- data.frame("logL" = rep(NA, n), df.model = rep(NA, n),
  	AIC = rep(NA, n), AICc = rep(NA, n), BIC = rep(NA, n),
  	Deviance = rep(NA, n), "Resid. df" = rep(NA, n),
  	Test = rep(NA, n), DevDiff = rep(NA, n), 
  	df = rep(NA, n), "P(> DevDiff)" = rep(NA, n), check.names = FALSE)
	
	mod <- vector(mode = "character", length = length(mods))
	
	# utility function to remove white spaces at the beginning and at the end of character strings
  tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)

  for(i in seq(n)){
  	#    fm <- get(nam[i])
    fm <- mods[[i]]
  	# model name and formula
  	# paste + tr are necessary when there are interaction terms in the formula
    fb <- paste(tr(deparse(fm$formula)), collapse = " ")
    ft <- deparse(fm$phi.formula)
    mod[i] <- paste(names(mods)[i], ": ", "Mu = ", fb, "; Phi = ", ft, sep = "")
  	# anova table
    df.model <- fm$df.model ; akic <- AIC(fm)
    dfres <- df.residual(fm)
    dfr$logL[i] <- fm$logL
    dfr$df.model[i]    <- df.model
    dfr$AIC[i]  <- akic[, 3]
    dfr$AICc[i] <- akic[, 4]
    dfr$BIC[i]  <- -2 * fm$logL + log(df.model + dfres) * df.model
    dfr[i, 6]   <- deviance(fm)
    dfr[i, 7]   <- dfres
  }
  
	dfr$Test <- c(NA, paste(names(mods)[-nrow(dfr)], names(mods)[-1], sep = "-"))
  dfr$DevDiff  <- c(NA, 2 * diff(dfr[, 1]))
  dfr$df <- c(NA, diff(dfr[,2]))
  
	for(i in 2:n)
    dfr[i , 11] <- 1 - pchisq(abs(dfr$DevDiff[i]), df = abs(dfr$df[i]))
  rownames(dfr) <- names(mods)
	
	structure(list(models = mod, anova.table = dfr), class = "anova.aodml")

}
