fitted.aodml <- function(object, ..., what = c("mu", "nu", "eta", "phi")) {

  what <- match.arg(what)
  dat <- object$dat
	mu.f <- object$formula
	phi.f <- object$phi.formula
	modmatrix.b <- object$modmatrix.b
	modmatrix.phi <- object$modmatrix.phi
  offset <- object$offset
	b <- object$b
	phi <- object$phi
	
	# fitted mu
	nu <- as.vector(modmatrix.b %*% b)
  eta <- if(is.null(offset)) nu else nu + offset
	mu <- invlink(eta, type = object$link)
	
	# fitted phi
	zphi <- as.vector(modmatrix.phi %*% phi)

	switch(what, mu = mu, nu = nu, eta = eta, phi = zphi)
	}

fitted.aodql <- function(object, ...) fitted(object$fm, ...)



