\name{drs}
\alias{drs}
% \alias{drs-class}
\alias{print.drs}
\encoding{latin1}

\title{Test of Proportion Homogeneity between Groups using Donner's and Rao-Scott's Adjustments}

\description{The function  tests the homogeneity of probabilities between \eqn{J} groups (H_0: \eqn{\mu_1 = \mu_2 = ... = \mu_J}) from clustered binomial data \{\eqn{(n_1, m_1), (n_2, m_2), ..., (n_N, m_N)}\}, where \eqn{n_i} is the size of cluster \eqn{i}, \eqn{m_i} the number of \dQuote{successes} (proportions are \eqn{y = m/n}), and \eqn{N} the number of clusters. The function uses adjusted chi-squared statistics, with either the correction proposed by proposed by Donner (1989) or the correction proposed by Rao and Scott (1993).}

\usage{
  drs(formula, data, type = c("d", "rs"), C = NULL, pooled = FALSE)
  
  \method{print}{drs}(x, ...)
  }

\arguments{
  \item{formula}{An formula where the left-hand side is a matrix of the form \code{cbind(m, n-m)} (the modelled proportion is \eqn{m / n}). The right-hand side must specify a single grouping variable.}

  \item{type}{A character string: either \dQuote{d} for the Donner's test and \dQuote{rs} for the Rao and Scott's test.}

  \item{data}{A data frame containing \code{n}, \code{m}) and the grouping variable.}

  \item{C}{An optional vector of a priori \eqn{J} cluster correction factors used for the Donner's test or design effects factors used for the Rao-Scott's test. If C is set no NULL (default), it is calculated internally (see details).}

  \item{pooled}{Logical indicating if a pooled design effect is estimated over the \eqn{J} groups for the  Rao-Scott's test (see details). Default to FALSE.}

  \item{x}{An object of class \dQuote{drf}.}
  
  \item{...}{Further arguments to be passed to \code{print}.}
  
  }

\details{
  \bold{Donner's test}\cr
  The chi-squared statistic is adjusted with the correction factor \eqn{C_j} computed in each group \eqn{j}. The test statistic is given by:
  
  \deqn{X^2 = \sum_{j} ( (m_j - n_j * \mu)^2 / (C_j * n_j * \mu * (1 - \mu)) )}
  
where \eqn{\mu = \sum_{j} (m_j) / \sum_{j} (n_j)} and \eqn{C_j = 1 + (n_{A,j} - 1) * \rho}. \eqn{n_{A,j}} is a scalar depending on the cluster sizes, and \eqn{\rho} is the ANOVA estimate of the intra-cluster correlation assumed common across groups (see Donner, 1989 or Donner et al., 1994). The statistic is compared to a chi-squared distribution with \eqn{J - 1} degrees of freedom. Fixed correction factors \eqn{C_j} can be specified with the argument \code{C}.\cr

  \bold{Rao ans Scott's test}\cr
  The method uses design effects and \dQuote{effective} sample sizes. The design effect \eqn{C_j} in each group \eqn{j} is estimated by \eqn{C_j = v_{ratio,j} / v_{bin,j}}, where \eqn{v_{ratio,j}} is the variance of the ratio estimate of the probability in group \eqn{i} (Cochran, 1999, p. 32 and p. 66) and \eqn{v_{bin,j}} is the standard binomial variance. The \eqn{C_j} are used  to compute the effective sample sizes \eqn{n_{adj,j} = n_j / C_j}, the effective numbers of successes \eqn{m_{adj,j} = m_j / C_j} in each group \eqn{j}, and the overall effective proportion \eqn{mu_adj = \sum_{j} m_{adj,j} / \sum_{j} C_j}. The test statistic is obtained by substituting these quantities in the usual chi-squared statistic, yielding:
  
  \deqn{X^2 = \sum_{j} ( (m_{adj,j} - n_{adj,j} * muadj)^2 / (n_{adj,j} * muadj * (1 - muadj)) )}

which is compared to a chi-squared distribution with \eqn{J - 1} degrees of freedom.\cr
A pooled design effect over the \eqn{J} groups is estimated if argument \code{pooled = TRUE} (see Rao and Scott, 1993, Eq. 6). Fixed design effects \eqn{C_j} can be specified with the argument \code{C}.
}

\value{An object of class \code{drs}, printed with \code{print.drs}.}

\references{
  Donner, A., 1989. \emph{Statistical methods in ophthalmology: an adjusted chi-squared approach}. Biometrics 45, 605-611.\cr
  Donner, A., 1993. \emph{The comparison of proportions in the presence of litter effects}. Prev. Vet. Med. 18, 17-26.\cr
  Donner, A., Eliasziw, M., Klar, N., 1994. \emph{A comparison of methods for testing homogeneity of proportions in teratologic studies}. Stat. Med. 13, 1253-1264.
  }

\seealso{\code{\link{chisq.test}}}

\examples{
data(dja)
# Donner
drs(formula = cbind(m, n - m) ~ group, data = dja, type = "d")
# Rao and Scott
drs(formula = cbind(m, n - m) ~ group, type = "rs", data = dja)
drs(formula = cbind(m, n - m) ~ group, type = "rs", data = dja, pooled = TRUE)
# standard chi2 test
drs(formula = cbind(m, n - m) ~ group, data = dja, type = "d", C = c(1:1))
drs(formula = cbind(m, n - m) ~ group, data = dja, type = "rs", C = c(1:1))
}

\keyword{datagen}
