\name{splitbin}
\alias{splitbin}
\encoding{latin1}

\title{Split Grouped Data Into Individual Data}

\description{The function splits grouped binomial data and optional covariates to individual binary data. Two types of grouped data are managed by \code{splitbin}:
     \itemize{
        \item Grouped data with weights;
        \item Grouped data of form \code{cbind(success, failure)}.
        }
When weights, successes or failures involve non-integer numbers, these numbers are rounded (using \code{round()}) before splitting.}

\usage{splitbin(formula, data, id = "idbin")}

\arguments{
  \item{formula}{A formula. The left-hand side describes the grouped data. The right-hand side describes the covariates. See examples for syntax.}
  \item{data}{A data frame where all the variables described in \code{formula} are found.}
  \item{id}{An optional character string naming the identifier (= grouping factor). Default to \dQuote{idbin}.}
  }

\value{A data frame built according to the formula and function used in the call.}

\examples{
# grouped data with weights
z <- data.frame(
    m = c(0, 1, 0, 1),
    f1 = c("A", "A", "B", "B"),
    f2 = c("C", "D", "C", "D"),
    n = c(4, 2, 1, 3)
    )
z
splitbin(formula = n ~ f1, data = z)$tab
splitbin(formula = n ~ f1 + f2 + m , data = z)$tab

# grouped data of form "cbind(success, failure)"
z <- data.frame(
    m = c(4, 1),
    n = c(5, 3),
    f1 = c("A", "B"),
    f2 = c("C", "D")
    )
z
splitbin(formula = cbind(m, n - m) ~ 1, data = z)$tab
splitbin(formula = cbind(m, n - m) ~ f1 + f2, data = z)$tab
splitbin(formula = cbind(m, n - m) ~ f1 + f2, data = z)$tab
}

\keyword{datagen}

