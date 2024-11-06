\name{iccbin}
\alias{iccbin}
% \alias{iccbin-class}
\alias{print.iccbin}
\encoding{latin1}

\title{Intra-Cluster Correlation for Clustered Binomial data}

\description{The function estimates the intraclass correlation \eqn{\rho} from clustered binomial data:

\{\eqn{(n_1, m_1), (n_2, m_2), ..., (n_N, m_N)}\},

where \eqn{n_i} is the size of cluster \eqn{i}, \eqn{m_i} the number of \dQuote{successes} (proportions are \eqn{y = m/n}), and \eqn{N} the number of clusters. The function uses a one-way random effect model. Three estimates, corresponding to methods referred to as \dQuote{A}, \dQuote{B} and \dQuote{C} in Goldstein et al. (2002), can be returned.}

\usage{
  iccbin(n, m, method = c("A", "B", "C"), nAGQ = 1, M = 1000)
  
  \method{print}{iccbin}(x, ...)
  }

\arguments{
  \item{n}{A vector of the sizes of the clusters.}
  \item{m}{A vector of the numbers of successes (proportions are \eqn{y = m / n}).}
  \item{method}{A character (\dQuote{A}, \dQuote{B} or \dQuote{C}) defining the calculation method. See Details.}
  \item{nAGQ}{Same as in function \code{glmer} of package \pkg{lme4}. Only for methods \dQuote{A} and \dQuote{B}. Default to 1.}
  \item{M}{Number of Monte Carlo (MC) replicates used in method \dQuote{B}. Default to 1000.}
  \item{x}{An object of class \dQuote{iccbin}.}
  \item{...}{Further arguments to ba passed to \dQuote{print}.}
  }

\details{
Before computations, the clustered data are split to binary \dQuote{0/1} observations \eqn{y_{ij}} (observation \eqn{j} in cluster \eqn{i}). The methods of calculation are described in Goldstein et al. (2002).\cr\cr

Methods "A" and "B" use the 1-way logistic binomial-Gaussian model

  \deqn{y_{ij} | \mu_{ij} \sim Bernoulli(\mu_{ij})}
  
  \deqn{logit(\mu_{ij}) = b_0 + u_i,}

where \eqn{b_0} is a constant and \eqn{u_i} a cluster random effect with \eqn{u_i \sim N(0, s^2_u)}. The ML estimate of the variance component \eqn{s^2_u} is calculated with the function \code{glmer} of package \pkg{lme4}. The intra-class correlation \eqn{\rho = Corr[y_{ij}, y_{ik}]} is then calculated from a first-order model linearization around \eqn{E[u_i]=0} in method \dQuote{A}, and with Monte Carlo simulations in method \dQuote{B}.\cr\cr

Method "C" provides the common ANOVA (moment) estimate of \eqn{\rho}. For details, see for instance Donner (1986), Searle et al. (1992) or Ukoumunne (2002).
}

\value{An object of class \code{iccbin}, printed with \code{print.iccbin}.}

\references{
  Donner A., 1986, \emph{A review of inference procedures for the intraclass correlation coefficient in the one-way random effects model}. International Statistical Review 54, 67-82.\cr
  Searle, S.R., Casella, G., McCulloch, C.E., 1992. \emph{Variance components}. Wiley, New York.\cr
  Ukoumunne, O. C., 2002. \emph{A comparison of confidence interval methods for the intraclass correlation coefficient in cluster randomized trials}. Statistics in Medicine 21, 3757-3774.\cr
  Golstein, H., Browne, H., Rasbash, J., 2002. \emph{Partitioning variation in multilevel models}.  Understanding Statistics 1(4), 223-231.
  }

\seealso{\code{\link[lme4:glmer]{glmer}}\cr}

\examples{
data(rats)
z <- rats[rats$group == "TREAT", ]
# A: glmm (model linearization)
iccbin(z$n, z$m, method = "A")
iccbin(z$n, z$m, method = "A", nAGQ = 10)
# B: glmm (Monte Carlo)
iccbin(z$n, z$m, method = "B")
iccbin(z$n, z$m, method = "B", nAGQ = 10, M = 1500)
# C: lmm (ANOVA moments)
iccbin(z$n, z$m, method = "C")

\dontrun{
  # Example of CI calculation with nonparametric bootstrap
  require(boot)
  foo <- function(X, ind) {
    n <- X$n[ind]
    m <- X$m[ind]
    iccbin(n = n, m = m, method = "C")$rho
    }
  res <- boot(data = z[, c("n", "m")], statistic = foo, R = 500, sim = "ordinary", stype = "i")
  res
  boot.ci(res, conf = 0.95, type = "basic")
  }
  
}

\keyword{datagen}
