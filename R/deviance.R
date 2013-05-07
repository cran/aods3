deviance.aodml <- function(object, ...) 2 * (sum(object$lmax) - sum(object$l))

deviance.aodql <- function(object, ...) deviance(object$fm, ...)

