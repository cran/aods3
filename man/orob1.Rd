\name{orob1}
\alias{orob1}
\docType{data}
\encoding{latin1}

\title{Germination Data}

\description{Data describing the germination for seed \emph{Orobanche cernua} cultivated in three dilutions of a bean root extract (Crowder, 1978, Table 1). The mean proportions of the three groups are 0.142, 0.872 and 0.842, and the overall mean is 0.614.}

\usage{data(orob1)}

\format{
  A data frame with 16 observations on the following 3 variables.
  \describe{
    \item{dilution}{a factor with 3 levels: \code{1/1}, \code{1/25} and \code{1/625}}.
    \item{n}{a numeric vector: the number of seeds exposed to germination.}
    \item{m}{a numeric vector: the number of seeds which actually germinated.}
    }
  }

\references{
  Crowder, M.J., 1978. \emph{Beta-binomial anova for proportions}. Appl. Statist. 27, 34-37.
  }

\keyword{datasets}
