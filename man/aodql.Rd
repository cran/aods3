\name{aodql}
\alias{aodql}
\alias{anova.aodql}
\alias{coef.aodql}
\alias{deviance.aodql}
\alias{df.residual.aodql}
\alias{fitted.aodql}
\alias{logLik.aodql}
\alias{predict.aodql}
\alias{print.aodql}
\alias{residuals.aodql}
\alias{summary.aodql}
\alias{vcov.aodql}
\encoding{latin1}

\title{QL/MM Estimation of Generalized Linear Models for Overdispersed Count Data}

\description{From clustered data, the function fits generalized linear models containing an over-dispersion parameter \eqn{\Phi} using quasi-likelihood estimating equations for the mean \eqn{\mu} and a moment estimator for \eqn{\Phi}.\cr

For binomial-type models, data have the form \{\eqn{(n_1, m_1), (n_2, m_2), ..., (n_N, m_N)}\}, where \eqn{n_i} is the size of cluster \eqn{i}, \eqn{m_i} the number of \dQuote{successes}, and \eqn{N} the number of clusters. The response is the proportion \eqn{y = m/n}. \cr

For Poisson-type models, data can be of two forms. When modeling \dQuote{simple counts}, data have the form \{\eqn{m_1, m_2, ..., m_N}\}, where \eqn{m_i} is the number of occurences of the event under study. When modeling rates (e.g. hazard rates), data have the same form as for the BB model, where \eqn{n_i} is the denominator of the rate for cluster \eqn{i} (considered as an \dQuote{offset}, i.e. a constant known value) and \eqn{m_i} the number of occurences of the event. For both forms of data, the response is the count \eqn{y = m}. }

\usage{
  aodql(formula,
        data,
        family = c("qbin", "qpois"),
        link = c("logit", "cloglog", "probit"),
        method = c("chisq", "dev"),
        phi = NULL,
        tol = 1e-5, ...)
        
  \method{anova}{aodql}(object, ...)
  \method{coef}{aodql}(object, ...)
  \method{deviance}{aodql}(object, ...)
  \method{df.residual}{aodql}(object, ...)
  \method{fitted}{aodql}(object, ...)
  \method{logLik}{aodql}(object, ...)
  \method{predict}{aodql}(object, ...)
  \method{print}{aodql}(x, ...)
  \method{residuals}{aodql}(object, ...)
  \method{summary}{aodql}(object, ...)
  \method{vcov}{aodql}(object, ...)  
  }

\arguments{
  
  \item{formula}{A formula for the mean \eqn{\mu}, defining the parameter vector \eqn{b} (see details). For binomial-type models, the left-hand side of the formula must be of the form \code{cbind(m, n - m) ~ ...} where the fitted proportion is \code{m/n}. For Poisson-type models, the left-hand side of the formula must be of the form \code{m ~ ...} where the fitted count is \code{m}. To fit a rate, argument \code{offset} must be used in the right-hand side of the formula (see examples).}
  
  \item{data}{A data frame containing the response (\code{m} and, optionnally, \code{n}) and the explanatory variable(s).}
    
  \item{family}{Define the model which is fitted: \dQuote{qbin} for binomial-type models and \dQuote{qpois} for Poisson-type models.}
  
  \item{link}{For binomial-type models only. Define the link function \eqn{g} for the mean \eqn{\mu}: \dQuote{cloglog}, \dQuote{logit} (default) or \dQuote{probit}. For Poisson-type models, \code{link} is automatically set to \dQuote{log}.}

  \item{method}{For function \code{aodql}, define the statistics used for the moment estimation of \eqn{phi}; legal values are \dQuote{chisq} (default) for the chi-squared statistics or \dQuote{dev} for the deviance statistics.}  
  
  \item{phi}{An optional value defining the over-dispersion parameter \eqn{\Phi} if it is set as constant. Default to NULL (in that case, \eqn{\Phi} is estimated).}

  \item{tol}{A positive scalar (default to 0.001). The algorithm stops at iteration \eqn{r + 1} when the condition
  \eqn{\chi{^2}[r+1] - \chi{^2}[r] <= tol} is met by the \eqn{\chi^2}{chi-squared} statistics .}
  
  \item{\dots}{Further arguments to passed to the appropriate functions.}

  \item{object}{An object of class \dQuote{aodql}.}

  \item{x}{An object of class \dQuote{aodql}.}
  }

\details{

\bold{Binomial-type models}\cr

For a given cluster \eqn{(n, m)}, the model is

  \deqn{m | \lambda,n \sim Binomial(n, \lambda)}

where \eqn{\lambda} follows a random variable of mean \eqn{E[\lambda] = \mu} and variance \eqn{Var[\lambda] = \Phi * \mu * (1 - \mu)}. The marginal mean and variance of \eqn{m} are

  \deqn{E[m] = n * \mu}

  \deqn{Var[m] = n * \mu * (1 - \mu) * (1 + (n - 1) * \Phi)}
    
The response in \code{aodql} is \eqn{y = m/n}.  The mean is \eqn{E[y] = \mu}, defined such as \eqn{\mu = g^{-1}(X * b) =  g^{-1}(\nu)}, where \eqn{g} is the link function, \eqn{X} is a design-matrix, \eqn{b} a vector of fixed effects and \eqn{\nu = X * b} is the corresponding linear predictor. The variance is \eqn{Var[y] = (1 / n) * \mu * (1 - \mu) * (1 + (n - 1) * \Phi)}.\cr 

\bold{Poisson-type models}\cr

\emph{------ Simple counts (model with no offset)}\cr

For a given cluster \eqn{(m)}, the model is

  \deqn{y | \lambda \sim Poisson(\lambda)}

where \eqn{\lambda} follows a random distribution of mean \eqn{\mu} and variance \eqn{\Phi * \mu^2}. The mean and variance of the marginal distribution of \eqn{m} are

  \deqn{E[m] = \mu}

  \deqn{Var[m] = \mu + \Phi * \mu^2}

The response in \code{aodql} is \eqn{y = m}. The mean is \eqn{E[y] = \mu}, defined such as \eqn{\mu = exp(X * b) =  exp(\nu)}. The variance is \eqn{Var[y] = \mu + \Phi * \mu^2}.\cr

\emph{------ Rates (model with offset)}\cr

For a given cluster \eqn{(n, m)}, the model is

  \deqn{m | \lambda,n \sim Poisson(\lambda)}

where \eqn{\lambda} follows the same random distribution as for the case with no offset. The marginal mean and variance are

  \deqn{E[m | n] = \mu}

  \deqn{Var[m | n] = \mu + \Phi * \mu^2}

The response in \code{aodql} is \eqn{y = m}. The mean is \eqn{E[y] = \mu}, defined such as \eqn{\mu = exp(X * b + log(n)) =  exp(\nu + log(n)) = exp(\eta)}, where \eqn{log(n)} is the offset. The variance is \eqn{Var[y] = \mu + \Phi * \mu^2}.\cr

\bold{Other details}\cr

Vector \eqn{b} and parameter \eqn{\Phi} are estimated iteratively, using procedures referred to as "Model I" in Williams (1982) for binomial-type models, and "Procedure II" in Breslow (1984) for Poisson-type models.\cr

Iterations are as follows. Quasi-likelihood estimating equations (McCullagh & Nelder, 1989) are used to estimate \eqn{b} (using function \code{glm} and its \code{weights} argument), \eqn{\Phi} being set to a constant value. Then, \eqn{\Phi} is calculated by the moment estimator, obtained by equalizing the goodness-of-fit statistic (chi-squared \code{X2} or deviance \code{D}) of the model to its degrees of freedom.\cr
Parameter \eqn{\Phi} can be set as constant, using argument \code{phi}. In that case, only \eqn{b}  is estimated.
}

\value{An object of class \code{aodql}, printed and summarized by various functions.}

\references{
  Breslow, N.E., 1984. \emph{Extra-Poisson variation in log-linear models}. Appl. Statist. 33, 38-44.\cr
  Moore, D.F., 1987, \emph{Modelling the extraneous variance in the presence of extra-binomial variation}. 
    Appl. Statist. 36, 8-14.\cr
  Moore, D.F., Tsiatis, A., 1991. \emph{Robust estimation of the variance in moment methods for extra-binomial
  and extra-poisson variation}. Biometrics 47, 383-401.
  McCullagh, P., Nelder, J. A., 1989, 2nd ed. \emph{Generalized linear models}. New York, USA: Chapman and Hall.\cr
  Williams, D.A., 1982, \emph{Extra-binomial variation in logistic linear models}. Appl. Statist. 31, 144-148.
  }
  
\seealso{\code{\link{glm}}}
  
\examples{

#------ Binomial-type models

data(orob2)
fm <- aodql(cbind(m, n - m) ~ seed, data = orob2, family = "qbin")
coef(fm)
vcov(fm)
summary(fm)
# chi2 tests of the seed factor in fm
wald.test(b = coef(fm), varb = vcov(fm), Terms = 2)

# chi-2 vs. deviance statistic to estimate phi
fm1 <- aodql(cbind(m, n - m) ~ seed + root, data = orob2, family = "qbin")
fm2 <- aodql(cbind(m, n - m) ~ seed + root, data = orob2, family = "qbin", method = "dev")
coef(fm1)
coef(fm2)
fm1$phi
fm2$phi
vcov(fm1)
vcov(fm2)
gof(fm1)
gof(fm2)

# estimate with fixed phi
fm <- aodql(cbind(m, n - m) ~ seed, data = orob2, family = "qbin", phi = 0.05)
coef(fm)
vcov(fm)
summary(fm)

#------ Poisson-type models

data(salmonella)
fm <- aodql(m ~ log(dose + 10) + dose, data = salmonella, family = "qpois")
coef(fm)
vcov(fm)
summary(fm)
# chi2 tests of the "log(dose + 10) + dose" factors
wald.test(b = coef(fm), varb = vcov(fm), Terms = 2:3)

# chi-2 vs. deviance statistic to estimate phi
fm1 <- aodql(m ~ log(dose + 10) + dose, data = salmonella, family = "qpois")
fm2 <- aodql(m ~ log(dose + 10) + dose, data = salmonella, family = "qpois", method = "dev")
coef(fm1)
coef(fm2)
fm1$phi
fm2$phi
vcov(fm1)
vcov(fm2)
gof(fm1)
gof(fm2)

# estimate with fixed phi
fm <- aodql(m ~ log(dose + 10) + dose, data = salmonella, family = "qpois", phi = 0.05)
coef(fm)
vcov(fm)
summary(fm)

# modelling a rate
data(dja)
# rate "m / trisk"
fm <- aodql(formula = m ~ group + offset(log(trisk)), data = dja, family = "qpois")
summary(fm)

}

\keyword{regression}
