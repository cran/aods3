predict.aodml <- function(object, ..., type = c("link", "response"), se.fit = FALSE, newdata = NULL) {
  type <- match.arg(type)
		
	dat <- object$dat
	mu.f <- object$formula
	phi.f <- object$phi.formula
	modmatrix.b <- object$modmatrix.b
	modmatrix.phi <- object$modmatrix.phi
  offset <- object$offset
		
	if(!is.null(newdata)) {
		modmatrix.b <- model.matrix(object = mu.f, data = newdata)
		mf <- model.frame(formula = mu.f, data = newdata)
  	offset <- model.offset(mf)
	}
	
	X <- modmatrix.b
	b <- object$b

	# predict eta ("link")
	nu <- as.vector(X %*% b)
  eta <- if(is.null(offset)) nu else nu + offset
  varb <- vcov(object)
  vareta <- X %*% varb %*% t(X)
  se.eta <- as.vector(sqrt(diag(vareta)))
	
  # predict mu ("response")
	mu <- invlink(eta, type = object$link)
  J <- switch(
  	object$link,
    cloglog = diag(-(1 - mu) * log(1 - mu), nrow = length(mu)),
    log = diag(mu, nrow = length(mu)),
    logit = diag(mu * (1 - mu), nrow = length(mu)),
  	probit = diag(dnorm(eta), nrow = length(mu))
  	)
  varmu <- J %*% vareta %*% J
  se.mu <- as.vector(sqrt(diag(varmu)))
  
	if(!se.fit)
    res <- switch(type, response = mu, link = eta)
  else
    res <- switch(type, response = list(fit = mu, se.fit = se.mu), link = list(fit = eta, se.fit = se.eta))
	res
	}

predict.aodql <- function(object, ...) predict(object$fm, ...)
