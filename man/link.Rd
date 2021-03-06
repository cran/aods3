\name{link}
\alias{link}
\encoding{latin1}

\title{Transformation from the Observation Scale to the Link Scale}

\description{The function transforms a variable from the observation scale (probability or count) to the link scale.}

\usage{link(x, type = c("cloglog", "log", "logit", "probit"))}

\arguments{
  \item{x}{A vector of real numbers.}
  \item{type}{A character string: \dQuote{cloglog}, \dQuote{log}, \dQuote{logit} or \dQuote{probit}.}
  }

\value{
  clog-log: \eqn{log(-log(1 - x ))}\cr
  log: \eqn{log(x)}\cr
  logit: \eqn{log(x / (1 - x))}\cr
  probit: \eqn{qnorm(x)}
  }

\seealso{\code{\link{invlink}}}

\examples{
x <- seq(.001, .999, length = 100)
plot(x, link(x, type = "logit"), type = "l", lwd = 2, ylab = "link(proba.)")
lines(x, link(x, type = "cloglog"), lty = 2, lwd = 2)
grid(col = "black")
legend(0, 6, legend = c("logit(x)", "cloglog(x)"), lty = c(1, 2), bg = "white")
}

\keyword{math}
