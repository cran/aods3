\name{rabbits}
\alias{rabbits}
\docType{data}
\encoding{latin1}

\title{Rabbits Foetuses Survival Experiment}

\description{Experimental data for analyzing the effect of an increasing dose of a compound on the proportion of live foetuses affected (Paul, 1982, Table 1). Four treatment-groups were considered: control \dQuote{C}, low dose \dQuote{L}, medium dose \dQuote{M} and high dose \dQuote{H}. The animal species used in the experiment was banded Dutch rabbit.}

\usage{data(rabbits)}

\format{
  A data frame with 84 observations on the following 3 variables.
  \describe{
    \item{group}{a factor with levels \code{C}, \code{H}, \code{L} and \code{M}}
    \item{n}{a numeric vector: the total number of foetuses.}
    \item{m}{a numeric vector: the number of affected foetuses.}
    }
  }

\references{
  Paul, S.R., 1982. \emph{Analysis of proportions of affected foetuses in teratological experiments.} 
  Biometrics 38, 361-370.
  }

\keyword{datasets}
