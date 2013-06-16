\name{aodml}
\alias{aodml}
\alias{print.aodml}
\alias{summary.aodml}
\alias{print.summary.aodml}
\alias{anova.aodml}
\alias{print.anova.aodml}
\alias{fitted.aodml}
\alias{residuals.aodml}
\alias{coef.aodml}
\alias{logLik.aodml}
\alias{deviance.aodml}
\alias{df.residual.aodml}
\alias{AIC.aodml}
\alias{vcov.aodml}
\alias{predict.aodml}

\encoding{latin1}

\title{ML Estimation of Generalized Linear Models for Overdispersed Count Data}

\description{The function fits a beta-binomial (BB) or a negative binomial (NB) generalized linear model from clustered data. \cr

For the BB model, data have the form \{\eqn{(n_1, m_1), (n_2, m_2), ..., (n_N, m_N)}\}, where \eqn{n_i} is the size of cluster \eqn{i}, \eqn{m_i} the number of \dQuote{successes}, and \eqn{N} the number of clusters. The response is the proportion \eqn{y = m/n.} 

For the NB model, data can be of two types. When modeling simple counts, data have the form \{\eqn{m_1, m_2, ..., m_N}\}, where \eqn{m_i} is the number of occurences of the event under study. When modeling rates (e.g. hazard rates), data have the same form as for the BB model, where \eqn{n_i} is the denominator of the rate for cluster \eqn{i} (considered as an \dQuote{offset}, i.e. a constant known value) and \eqn{m_i} the number of occurences of the event. For both types of data, the response is the count \eqn{y = m}. }

\usage{
  aodml(formula,
  data,
  family = c("bb", "nb"),
  link = c("logit", "cloglog", "probit"), 
	phi.formula = ~ 1,
  phi.scale = c("identity", "log", "inverse"),
  phi.start = NULL,
  fixpar = list(),
  hessian = TRUE,
  method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
  control = list(maxit = 3000, trace = 0), ...)
  
  \method{print}{aodml}(x, ...)
  \method{summary}{aodml}(object, ...)
  \method{anova}{aodml}(object, ...)
  \method{print}{anova.aodml}(x, digits, ...)
  \method{fitted}{aodml}(object, ..., what = c("mu", "nu", "eta", "phi"))
  \method{residuals}{aodml}(object, ..., type = c("deviance", "pearson", "response"))
  \method{coef}{aodml}(object, ...)
  \method{df.residual}{aodml}(object, ...)
  \method{logLik}{aodml}(object, ...)
  \method{deviance}{aodml}(object, ...)
  \method{AIC}{aodml}(object, ..., k = 2)
  \method{vcov}{aodml}(object, ...)
  \method{predict}{aodml}(object, ..., type = c("link", "response"), se.fit = FALSE, newdata = NULL)
  }

\arguments{

  \item{formula}{A formula for the mean \eqn{\mu}, defining the parameter vector \eqn{b} (see details).
  
For the BB model, the left-hand side of the formula must be of the form 

\code{cbind(m, n - m) ~ ...}

where the fitted proportion is \code{m/n}.

For the NB model, the left-hand side of the formula must be of the form 

\code{m ~ ...}

where the fitted count is \code{m}. To fit a rate, argument \code{offset} must be used in the right-hand side of the formula (see examples).}
  
  \item{data}{A data frame containing the response (\code{m} and, optionnally, \code{n}) and the explanatory variable(s).}
  
  \item{family}{Define the model which is fitted: \dQuote{bb} for the BB model and \dQuote{nb} for the NB model.}
  
  \item{link}{For the BB model only. Define the link function \eqn{g} for the mean \eqn{\mu}: \dQuote{cloglog}, \dQuote{logit} (default) or \dQuote{probit}. For the NB model, \code{link} is automatically set to \dQuote{log}.}
  
  \item{phi.formula}{A right-hand side formula to model optional heterogeneity for the over-dispersion parameter \eqn{\Phi} (see details). Only one single factor is allowed.
  
Default to \code{formula(~ 1)} (i.e. no heterogeneity).}
  
  \item{phi.scale}{Scale on which \eqn{\Phi} is estimated (see details): \dQuote{identity} (default), \dQuote{log} or \dQuote{inverse}.}

  \item{phi.start}{Optional starting values for \eqn{\Phi}. Default to 0.1.}
  
  \item{fixpar}{An optional list of 2 vectors of same length (\eqn{>=1}) to set some parameters as constant in the model.
  
The first vector indicates which parameters are set as constant (i.e., not optimized) in the global parameter vector \eqn{(b, \Phi)}.

The second vector indicates the corresponding values. For instance,

\code{fixpar = list(c(4, 5), c(0, 0.1))}

means that the 4th and 5th components of vector \eqn{(b, \Phi)} are set to 0 and 0.1. Argument \code{fixpar} can be usefull, for instance, to calculate profiled log-likehoods.}

  \item{hessian}{A logical (default to \code{TRUE}). If \code{FALSE}, the hessian and the variances-covariances matrices of the parameters are not calculated.}

  \item{method}{Define the method used for the optimization (see \code{\link{optim}}).}
  
  \item{control}{A list to control the optimization parameters. See \code{\link{optim}}. By default, the maximum number of iterations is set to 3000, and trace is set to 0 to avoid spurious warnings.}
  
  \item{object}{An object of class \code{aodml}}

  \item{x}{An object of class \code{aodml}}
  
  \item{digits}{Number of digits to print in print.summary.aodml and print.anova.aodml.
                Default to \code{max(3, getOption("digits") - 3)}}

  \item{\dots}{Further arguments passed to \code{\link{optim}} (e.g. argument \code{method} if using function \code{aodml}), or further objects of class \code{aodml} (function \code{anova.aodml}), or further arguments passed to \code{print.aodml} and \code{print.anova.aodml}.}
  
  \item{what}{For function \code{fitted}, a character string indicating the type of fitted values to be returned: legal values are \dQuote{mu} for the fitted response; \dQuote{nu} for the fitted linear predictor without offset (link scale); \dQuote{eta} for the fitted linear predictor with offset (link scale); \dQuote{phi} for the fitted overdispersion coefficient.}

  \item{type}{For function \code{residuals}, a character string indicating the type of residuals to be computed; legal values are \dQuote{deviance} for the deviance's residuals, \dQuote{pearson} for the Pearson's residuals, and \dQuote{response} for the response. For function \code{predict}, a character string indicating the type of prediction to be computed; legal values are \dQuote{link} and \dQuote{response}.}
  
  \item{k}{Numeric scalar for the penalty parameter used to compute the information criterion. The default value (\eqn{k = 2}) is the regular AIC = -2 * logLik + 2 * \eqn{p}, where \eqn{p} is the number of model coefficients. NB: for AICc, \eqn{k} is set to 2, and AICc = AIC + 2 * \eqn{p} * \eqn{(p + 1)} / \eqn{(n - p - 1)}, with \eqn{n} the number of observations.}
  
  \item{se.fit}{Logical scalar indicating whether standard errors should be computed for the predicted values. Default to \code{FALSE}.}
  
  \item{newdata}{A \code{data.frame} containing the explanatory variables - and possibly the offset - for the values of which predictions are to be made.}
  }


\details{

\bold{Beta-binomial model (BB)}\cr

For a given cluster \eqn{(n, m)}, the model is

  \deqn{m | \lambda,n ~ Binomial(n, \lambda)}

where \eqn{\lambda} follows a Beta distribution \eqn{Beta(a_1, a_2)}. Noting \eqn{B} the beta function, the marginal (beta-binomial) distribution of \eqn{m} is

  \deqn{P(m | n) = C(n, m) * B(a_1 + m, a_2 + n - m) / B(a_1, a_2)}

Using the re-parameterization \eqn{\mu = a_1 / (a_1 + a_2)} and \eqn{\Phi = 1 / (a_1 + a_2 + 1)}, the marginal mean and variance are

  \deqn{E[m] = n * \mu}

  \deqn{Var[m] = n * \mu * (1 - \mu) * (1 + (n - 1) * \Phi)}

The response in \code{aodml} is \eqn{y = m/n}. The mean is \eqn{E[y] = \mu}, defined such as \eqn{\mu = g^{-1}(X * b) =  g^{-1}(\nu)}, where \eqn{g} is the link function, \eqn{X} is a design-matrix, \eqn{b} a vector of fixed effects and \eqn{\nu = X * b} is the corresponding linear predictor. The variance is \eqn{Var[y] = (1 / n) * \mu * (1 - \mu) * (1 + (n - 1) * \Phi)}.\cr

\bold{Negative binomial model (NB)}\cr

\emph{------ Simple counts (model with no offset)}\cr

For a given cluster \eqn{(m)}, the model is

  \deqn{y | \lambda ~ Poisson(\lambda)}

where \eqn{\lambda} follows a Gamma distribution of mean \eqn{\mu} and shape \eqn{k} (\eqn{Var[\lambda] = \mu^2 / k}). Noting \eqn{G} the gamma function, the marginal (negative binomial) distribution of \eqn{m} is

  \deqn{P(m) = {G(m+k) / (m! * G(k))} * (k / (k + \mu))^k * (\mu / (k + \mu))^m}

Using the re-parameterization \eqn{\Phi = 1 / k}, the marginal mean and variance are

  \deqn{E[m] = \mu}

  \deqn{Var[m] = \mu + \Phi * \mu^2}

The response in \code{aodml} is \eqn{y = m}. The mean is \eqn{E[y] = \mu}, defined such as \eqn{\mu = exp(X * b) =  exp(\nu)}. The variance is \eqn{Var[y] = \mu + \Phi * \mu^2}.\cr

\emph{------ Rates (model with offset)}\cr

For a given cluster \eqn{(n, m)}, the model is

  \deqn{m | \lambda,n ~ Poisson(\lambda)}

The marginal (negative binomial) distribution \eqn{P(m|n)} is the same as for the case with no offset (\eqn{= P(m)}). The response in \code{aodml} is \eqn{y = m}. The mean is \eqn{E[y] = \mu}, defined such as \eqn{\mu = exp(X * b + log(n)) =  exp(\nu + log(n)) = exp(\eta)}, where \eqn{log(n)} is the offset. The variance is \eqn{Var[y] = \mu + \Phi * \mu^2}.\cr

\bold{Other details}\cr

Argument \code{phi.scale} of function \code{aodml} enables to estimate the over-dispersion parameter under different scales.

If \code{phi.scale = "identity"} (Default), the function estimates \eqn{\Phi}. 

If \code{phi.scale = "log"}, the function estimates \eqn{log(\Phi)}.

If \code{phi.scale = "inverse"}, the function estimates \eqn{1/\Phi}.\cr

The full parameter vector returned by \code{aodml}, say \code{param}, is equal to \eqn{(b, \Phi)}. This vector is estimated by maximizing the log-likelihood of the marginal model using function \code{optim}. The estimated variances-covariances matrix of \code{param} is calculated by the inverse of the observed hessian matrix returned by \code{optim}, and is referred to as \code{varparam}.

}

\value{An object of class \code{aodml}, printed and summarized by various functions. Function \code{deviance.aodml} returns the value \code{-2 * (logL - logL_max)}. The \dQuote{deviance} used in function \code{AIC.aodml} to calculate AIC and AICc is \code{-2 * logL}.}

\references{
  Crowder, M.J., 1978. \emph{Beta-binomial anova for proportions}. Appl. Statist. 27, 34-37.\cr
  Griffiths, D.A., 1973. \emph{Maximum likelihood estimation for the beta-binomial distribution and an application
  to the household distribution of the total number of cases of disease}. Biometrics 29, 637-648.\cr
  Lawless, J.F., 1987. \emph{Negative binomial and mixed Poisson regression}. The Canadian Journal of Statistics, 15(3): 209-225.\cr
  McCullagh, P., Nelder, J. A., 1989, 2nd ed. \emph{Generalized linear models}. New York, USA: Chapman and Hall.\cr
  Prentice, R.L., 1986. \emph{Binary regression using an extended beta-binomial distribution, with discussion of
  correlation induced by covariate measurement errors}. J.A.S.A. 81, 321-327.\cr
  Williams, D.A., 1975. \emph{The analysis of binary responses from toxicological experiments involving
  reproduction and teratogenicity}. Biometrics 31, 949-952.
  }
  
\seealso{\code{\link{glm}} and \code{\link{optim}}}
  
\examples{

#------ Beta-binomial model

data(orob2)
fm1 <- aodml(cbind(m, n - m) ~ seed, data = orob2, family = "bb")

# summaries
fm1
summary(fm1)
coef(fm1)
vcov(fm1)
logLik(fm1)
deviance(fm1)
AIC(fm1)
gof(fm1)

# predictions
cbind(
  fitted(fm1),
  fitted(fm1, what = "nu"),
  fitted(fm1, what = "eta"),
  fitted(fm1, what = "phi")
)
predict(fm1, type = "response", se.fit = TRUE)
newdat <- data.frame(seed = c("O73", "O75"))
predict(fm1, type = "response", se.fit = TRUE, newdata = newdat)

# model with heterogeneity in phi
fm <- aodml(cbind(m, n - m) ~ seed, data = orob2,
  family = "bb", phi.formula = ~ seed)
summary(fm)
AIC(fm1, fm)

# various phi scales
fm <- aodml(cbind(m, n - m) ~ seed, data = orob2, family = "bb")
fm$phi
fm$phi.scale
fm <- aodml(cbind(m, n - m) ~ seed, data = orob2, family = "bb",
            phi.scale = "log")
fm$phi
fm$phi.scale
fm <- aodml(cbind(m, n - m) ~ seed, data = orob2, family = "bb",
            phi.scale = "inverse")
fm$phi
fm$phi.scale

### Models with coefficient(s) set as constant

# model with 1 phi coefficient, set as constant "0.02"
fm <- aodml(formula = cbind(m, n - m) ~ seed * root, data = orob2,
  family = "bb", fixpar = list(5, 0.02))
fm$param
fm$varparam

# model with 2 phi coefficients, with the first set as constant ~ "0"
fm <- aodml(formula = cbind(m, n - m) ~ seed * root, data = orob2,
  family = "bb", phi.formula = ~ seed, fixpar = list(5, 1e-15))
fm$param
fm$varparam

# model with 2 phi coefficients, with the first set as constant ~ "0",
# and the mu intercept (1rst coef of vector b) set as as constant "-0.5"
fm <- aodml(formula = cbind(m, n - m) ~ seed * root, data = orob2,
            family = "bb", phi.formula = ~ seed,
            fixpar = list(c(1, 5), c(-0.5, 1e-15)))
fm$param
fm$varparam
  
### Model tests

# LR tests - non-constant phi
fm0 <- aodml(cbind(m, n - m) ~ 1, data = orob2, family = "bb")
fm2 <- aodml(cbind(m, n - m) ~ seed + root, data = orob2, family = "bb")
fm3 <- aodml(cbind(m, n - m) ~ seed * root, data = orob2, family = "bb")
anova(fm0, fm1, fm2, fm3)

# LR tests - constant phi
# phi is assumed to be estimated from fm3
fm2.bis <- aodml(cbind(m, n - m) ~ seed  + root, data = orob2,
                 family = "bb", fixpar = list(4, fm3$phi))
LRstat <- 2 * (logLik(fm3) - logLik(fm2.bis))  
pchisq(LRstat, df = 1, lower.tail = FALSE)  
  
# Wald test of the seed factor in fm1
wald.test(b = coef(fm3), varb = vcov(fm3), Terms = 4)

#------ Negative binomial model

### Modelling counts

data(salmonella)
fm1 <- aodml(m ~ log(dose + 10) + dose, data = salmonella, family = "nb")
## fm1 <- aodml(m ~ log(dose + 10) + dose, data = salmonella, family = "nb",
##              method = "Nelder-Mead")

# summaries
fm1
summary(fm1)
coef(fm1)
vcov(fm1)
logLik(fm1)
deviance(fm1)
AIC(fm1)
gof(fm1)

# predictions
cbind(
  fitted(fm1),
  fitted(fm1, what = "nu"),
  fitted(fm1, what = "eta"),
  fitted(fm1, what = "phi")
)
predict(fm1, type = "response", se.fit = TRUE)
newdat <- data.frame(dose = c(20, 40))
predict(fm1, type = "response", se.fit = TRUE, newdata = newdat)

# various phi scales
fm <- aodml(m ~ log(dose + 10) + dose, data = salmonella, family = "nb")
fm$phi
fm$phi.scale
fm <- aodml(m ~ log(dose + 10) + dose, data = salmonella,
            family = "nb", phi.scale = "log")
fm$phi
fm$phi.scale
fm <- aodml(m ~ log(dose + 10) + dose, data = salmonella,
            family = "nb", phi.scale = "inverse")
fm$phi
fm$phi.scale
  
# LR and Wald tests of the "log(dose + 10) + dose" factors
fm0 <- aodml(m ~ 1, data = salmonella, family = "nb")
anova(fm0, fm1)
fm0.bis <- aodml(m ~ 1, data = salmonella, family = "nb",
                 fixpar = list(2, fm1$phi))
LRstat <- 2 * (logLik(fm1) - logLik(fm0.bis))  
pchisq(LRstat, df = 2, lower.tail = FALSE)  
wald.test(b = coef(fm1), varb = vcov(fm1), Terms = 2:3)

### Modelling a rate
  
data(dja)
# rate "m / trisk"
fm <- aodml(formula = m ~ group + offset(log(trisk)),
            data = dja, family = "nb")
summary(fm)

fm <- aodml(formula = m ~ group + offset(log(trisk)),
            phi.formula = ~ group, data = dja, family = "nb",
            phi.scale = "log")
summary(fm)
  
# model with 1 phi coefficient, set as constant "0.8"
fm <- aodml(formula = m ~ group + offset(log(trisk)), data = dja,
            family = "nb", phi.formula = ~1, fixpar = list(3, 0.8))
fm$param
fm$varparam

# model with 2 phi coefficients, with the first set as constant ~ "0" in the identity scale
fm <- aodml(formula = m ~ group + offset(log(trisk)), data = dja,
            family = "nb", phi.formula = ~ group, phi.scale = "log",
            fixpar = list(4, -15))
fm$param
fm$varparam

# model with 2 phi coefficients, with the first set as constant "0" in the identity scale,
# and the mu intercept coefficient (1rst coef of vector b) set as as constant "-0.5"
fm <- aodml(formula = m ~ group + offset(log(trisk)), data = dja,
  family = "nb", phi.formula = ~ group, phi.scale = "log",
  fixpar = list(c(1, 4), c(-0.5, -15)))
fm$param
fm$varparam
}

\keyword{regression}
