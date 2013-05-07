aodml <- function(formula, data,
	family = c("bb", "nb"), link = c("logit", "cloglog", "probit"), 
	phi.formula = ~ 1, phi.scale = c("identity", "log", "exp", "inverse"), phi.ini = NULL,
	fixpar = list(), hessian = TRUE, control = list(maxit = 3000), ...) {
	
	call <- match.call(expand.dots = FALSE)
	fam <- match.arg(family)
	phi.scale <- match.arg(phi.scale)

	mu.f <- formula
	phi.f <- phi.formula
  if(phi.f != ~ 1) phi.f <- update.formula(phi.f, ~ . - 1)

	dat <- data
	
	mf <- model.frame(formula = mu.f, data = dat)
	resp <- model.response(mf)
	if(fam == "bb") {
		fam0 <- eval(parse(text = paste("binomial(link =", link,")")))
		link <- match.arg(link)
		m <- resp[, 1]
		n <- rowSums(resp)
		y <- m / n
	}
	if(fam == "nb") {
		fam0 <- "poisson"
		link <- "log"
		y <- as.vector(resp)
	}
  offset <- model.offset(mf)
	
	fam0 <- ifelse(fam == "nb", "poisson", "binomial")
  fm <- glm(formula = mu.f, family = fam0, data = dat)
  
	# model matrices
  modmatrix.b <- model.matrix(object = mu.f, data = dat)
  modmatrix.phi <- model.matrix(object = phi.f, data = dat)
  nbb <- ncol(modmatrix.b)
	nbphi <- ncol(modmatrix.phi)
	# total number of parameters (including the eventual ones that are set as constant)
	nbpar <- nbb + nbphi
	
	# Parameters set as constant (fixpar)
	# Argument fixpar must be a list with two components:
	# the rank(s) and the value(s) of the parameters that are set as constant 
	idpar <- idest <- 1:nbpar
	fixp <- FALSE
	idfix <- valfix <- NULL
	if(!is.null(unlist(fixpar))) {
		fixp <- TRUE
		idfix <- fixpar[[1]]
		valfix <- fixpar[[2]]
	  idest <- idpar[-idfix]
	}

	# initial parameter values
  b <- fm$coefficients
	if(is.null(phi.ini)) {
		z <- 0.1
		val <- switch(phi.scale, "identity" = z, "log" = log(z), "exp" = exp(z), "inverse" = 1 / z)
	} else
		val <- phi.ini
  phi.ini <- rep(val, nbphi)
	param.ini <- c(b, phi.ini)

	# -logL
  
  loglik.bb.k <- function(m, n, mu, k) {
    a <- (k - 1) * mu
    b <- (k - 1) * (1 - mu)
    lbeta(a + m, b + n - m) - lbeta(a, b) + lchoose(n, m)
  }
  
  loglik.nb.k <- function(y, mu, k)
    y * log(mu / (mu + k)) + k * log(k / (mu + k)) + lgamma(k + y) - lgamma(k) - lfactorial(y)
  
  minuslogL <- function(param0){
  	
  	param <- vector(length = nbpar)
  	param[idfix] <- valfix
  	param[idest] <- param0
  	
  	b <- param[1:nbb]
    nu <- as.vector(modmatrix.b %*% b)
    mu <- if(is.null(offset)) invlink(nu, type = link) else invlink(nu + offset, type = link)
    zphi <- as.vector(modmatrix.phi %*% param[(nbb + 1):(nbb + nbphi)])
    
  	# If some components of "zphi" are set
  	# to 0 (if "identity"), -Inf (if "log"), 1 (if "exp") or Inf (if "inverse") in fixpar,
  	# the distribution is set to Poisson and k (=1/phi) is not used
    
    k <- switch(phi.scale, identity = 1 / zphi, log = 1 / exp(zphi), exp = 1 / log(zphi), inverse = zphi)
    cnd <- k == Inf

    l <- vector(length = length(mu))
    
    if(fam == "bb") {
      l[cnd] <- dbinom(x = m[cnd], size = n[cnd], prob = mu[cnd], log = TRUE)
      l[!cnd] <- loglik.bb.k(m[!cnd], n[!cnd], mu[!cnd], k[!cnd])
    }
    
    if(fam == "nb") {
      l[cnd] <- dpois(x = y[cnd], lambda = mu[cnd], log = TRUE)
      l[!cnd] <- loglik.nb.k(y[!cnd], mu[!cnd], k[!cnd])
  	}
      
    fn <- sum(l)
    if(!is.finite(fn)) fn <- -1e20
    
    -fn
    
    }

	# fit
  res <- optim(par = param.ini[idest], fn = minuslogL, hessian = hessian, control = control, ...)
  
	## Results
  
	#param
	param0 <- res$par
  param <- vector(length = nbpar)
  param[idfix] <- valfix
  param[idest] <- param0
  namb <- colnames(modmatrix.b)
  namphi <- paste("phi", colnames(modmatrix.phi), sep = ".")
  names(param) <- c(namb, namphi)
  
  b <- param[seq(along = namb)]
  phi <- param[-seq(along = namb)]

  # varparam
  is.singular <- function(X) qr(X)$rank < nrow(as.matrix(X))
	varparam <- matrix(NA, ncol = nbpar, nrow = nbpar)
  singular.H0 <- NA
	if(hessian) {
    H0 <- res$hessian
    singular.H0 <- is.singular(H0)
    if(!singular.H0)
			varparam[idest, idest] <- qr.solve(H0)
		else
			warning("The hessian matrix was singular.\n")
	}
  dimnames(varparam)[[1]] <- dimnames(varparam)[[2]] <- names(param)
  
  # log-likelihood contributions
  nu <- as.vector(modmatrix.b %*% b)
  mu <- if(is.null(offset)) invlink(nu, type = link) else invlink(nu + offset, type = link)
  zphi <- as.vector(modmatrix.phi %*% phi)
  k <- switch(phi.scale, identity = 1 / zphi, log = 1 / exp(zphi), exp = 1 / log(zphi), inverse = zphi)
  cnd <- k == Inf
  l <- lmax <- vector(length = length(mu))

  if(fam == "bb") {
    l[cnd] <- dbinom(x = m[cnd], size = n[cnd], prob = mu[cnd], log = TRUE)
    l[!cnd] <- loglik.bb.k(m[!cnd], n[!cnd], mu[!cnd], k[!cnd])
    lmax[cnd] <- dbinom(x = m[cnd], size = n[cnd], prob = y[cnd], log = TRUE)
    lmax[!cnd] <- loglik.bb.k(m[!cnd], n[!cnd], y[!cnd], k[!cnd])
  }
  
  if(fam == "nb") {
    l[cnd] <- dpois(x = y[cnd], lambda = mu[cnd], log = TRUE)
    l[!cnd] <- loglik.nb.k(y[!cnd], mu[!cnd], k[!cnd])
    lmax[cnd] <- dpois(x = y[cnd], lambda = y[cnd], log = TRUE)
    lmax[!cnd] <- loglik.nb.k(y[!cnd], y[!cnd], k[!cnd])
    }  
  
  l[is.na(l)] <- 0
  lmax[is.na(lmax)] <- 0  
	
  # other results
	# if fixpar is not null, df.model is lower than nbpar
	df.model <- length(param0)
  logL <- -res$value
  df.residual <- length(y) - df.model
  iterations <- res$counts[1]
  code <- res$convergence
  msg <- if(!is.null(res$message)) res$message else character(0)
	
  #dev <- -2 * (logL - logLmax)
	#aic <- -2 * logL + 2 * df.model 
	
  structure(
  	list(
  		call = call, family = fam, link = link, dat = dat,
  		formula = mu.f, phi.formula = phi.f, phi.scale = phi.scale,
  		modmatrix.b = modmatrix.b, modmatrix.phi = modmatrix.phi,
  		resp = resp, offset = offset,
  		param = param, b = b, phi = phi,
    	varparam = varparam,
  		nbpar = nbpar, df.model = df.model, df.residual = df.residual,
  		logL = logL, l = l, lmax = lmax,
  		iterations = iterations, code = code, msg = msg,
      singular.hessian = singular.H0
  		),
  	class = "aodml")
}
