\name{wald.test}
\alias{wald.test}
\alias{print.wald.test}
\encoding{latin1}

\title{Wald Test for Model Coefficients}

\description{The function returns a Wald chi-squared test or a \eqn{F} test for a vector of model coefficients (possibly of length one), given its variance-covariance matrix.}

\usage{
  wald.test(b, varb, Terms = NULL, L = NULL, H0 = NULL,  df = NULL, verbose = FALSE, ...)
  
  \method{print}{wald.test}(x, ..., digits = max(3, getOption("digits") - 3))
  }

\arguments{
  \item{b}{
    A vector of coefficients with their var-cov matrix \code{varb}. Coefficients \code{b} and var-cov matrix are usually extracted using appropriate \code{coef} and \code{vcov} functions.
    }
  \item{varb}{
    A var-cov matrix of coefficients \code{b} (see above).
    }
  \item{Terms}{
    An optional integer vector specifying which coefficients should be \emph{jointly} tested, using a Wald chi-squared test or a\eqn{F} test. The elements of \code{varb} correspond to the columns or rows of the var-cov matrix given in \code{varb}. Default is \code{NULL}.
    }
  \item{L}{
    An optional matrix conformable to \code{b}, such as its product with \code{b} i.e., \code{L \%*\% b} gives the linear combinations of the coefficients to be tested. Default is \code{NULL}.
    }
  \item{H0}{
    A numeric vector giving the null hypothesis \eqn{H_0} for the test. It must be as long as \code{Terms} or must have the same number of columns as \code{L}. Default to 0 for all the coefficients to be tested.
    }
  \item{df}{
    A numeric vector giving the degrees of freedom to be used in an \eqn{F} test, i.e. the degrees of freedom of the residuals of the model from which \code{b} and \code{varb} were fitted. Default to NULL, for no \eqn{F} test. See the section \bold{Details} for more information.
    }
  \item{verbose}{
    A logical scalar controlling the amount of output information. The default is \code{FALSE}, providing minimum output.
    }
  \item{x}{An object of class \dQuote{wald.test}}
  \item{digits}{A numeric scalar indicating the number of digits to be kept after the decimal place.}
  \item{...}{
    Additional arguments to \code{print}.
    }
  }

\details{
  The assumption is that the coefficients follow asymptotically a multivariate normal distribution with mean equal to the model coefficients \code{b} and variance equal to their var-cov matrix \code{varb}.\cr\cr
  
  One (and only one) of \code{Terms} or \code{L} must be given. When \code{L} is given, it must have the same number of columns as the length of \code{b}, and the same number of rows as the number of linear combinations of coefficients.\cr\cr
  
  When \code{df} is given, the chi-squared Wald statistic is divided by \code{m}, the number of linear combinations of coefficients to be tested (i.e., \code{length(Terms)} or \code{nrow(L)}). Under the null hypothesis \eqn{H_0}, this new statistic follows an \eqn{F(m, df)} distribution.
  }

\value{An object of class \code{wald.test}, printed with \code{print.wald.test}.}

\references{
  Diggle, P.J., Liang, K.-Y., Zeger, S.L., 1994. \emph{Analysis of longitudinal data}. Oxford, Clarendon Press, 253 p.\cr
  Draper, N.R., Smith, H., 1998. \emph{Applied Regression Analysis}. New York, John Wiley & Sons, Inc., 706 p.
  }

\examples{
data(orob2)
fm <- aodql(cbind(m, n - m) ~ seed * root, data = orob2, family = "qbin")
# Wald chi2 test for the effect of root
wald.test(b = coef(fm), varb = vcov(fm), Terms = 3:4)
L <- matrix(c(0, 0, 1, 0, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
wald.test(b = coef(fm), varb = vcov(fm), L = L)
}

\keyword{htest}
