\name{invlink}
\alias{invlink}
\encoding{latin1}

\title{Transformation from the Link Scale to the Observation Scale}

\description{The function transforms a variable from the link scale to the observation scale (probability or count).}

\usage{invlink(x, type = c("cloglog", "log", "logit", "probit"))}

\arguments{
  \item{x}{A vector of real numbers.}
  \item{type}{A character string: \dQuote{cloglog}, \dQuote{log}, \dQuote{logit} or \dQuote{probit}.}
  }

\value{
  clog-log: \eqn{1 - exp(-exp(x))}\cr
  log: \eqn{exp(x)}\cr
  logit: \eqn{exp(x) / (1 + exp(x))}\cr
  probit: \eqn{pnorm(x)}
  }

\seealso{\code{\link{link}}}

\examples{
x <- seq(-5, 5, length = 100)
plot(x, invlink(x, type = "logit"), type = "l", lwd = 2, ylab = "Probability")
lines(x, invlink(x, type = "cloglog"), lty = 2, lwd = 2)
grid(col = "black")
legend(-5, 1, legend = c("alogit(x)", "acloglog(x)"), lty = c(1, 2), bg = "white")
}

\keyword{math}
