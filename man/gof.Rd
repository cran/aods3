\name{gof}
\alias{gof}
\alias{gof.default}
%\alias{gof-class}
%\alias{gof.default}
\alias{print.gof}
\encoding{latin1}

\title{Test of Goodness-of-Fit of Models for Count data}

\description{The function returns a chi-squared test of goodness of fit for models of class \code{glm}, \code{aodml} or \code{aodql}.}

\usage{
  gof(object)
  gof.default(object)
  
  \method{print}{gof}(x, ..., digits =  max(3, getOption("digits") - 3))
  }

\arguments{
  \item{object}{An object of class \code{glm}, \code{aodml} or \code{aodquasi}.}
  \item{x}{An object of class \code{gof}.}
  \item{digits}{A numerical scalar indicating the number of digits to be printed after the decimal place.}
  \item{...}{Further arguments passed to \code{print}.}
  }

\details{
Function \code{gof} calculates the deviance \eqn{D} and the Pearson chi-squared \eqn{X^2} statistics for the model under consideration. Let \eqn{y} be the observed response, and \eqn{E[y] = \mu} and \eqn{Var[y]} its mean and variance estimated from the model, statistic \eqn{X^2} is calculated by:

  \deqn{X^2 = \sum_{i}( (y_i - \mu)^2/Var[y_i] )}

Assuming that the data length is \eqn{N} and the number of the parameters in the model is \eqn{p}, eqn{D} and eqn{X^2} are compared to a chi-squared distribution with \eqn{N-p} degrees of freedom.

}

\value{An object of class \code{gof}, printed with \code{print.gof}.}

\references{Agresti, A. \emph{Categorical data analysis}. Wiley, 1990.}

\seealso{\code{\link{residuals}}, \code{\link{chisq.test}}}

\examples{
data(orob2)
fm1 <- glm(cbind(m, n - m) ~ seed, data = orob2, family = binomial)
fm2 <- aodml(cbind(m, n - m) ~ seed, data = orob2, family = "bb")
gof(fm1)
gof(fm2)
}

\keyword{datagen}
