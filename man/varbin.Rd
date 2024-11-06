\name{varbin}
\alias{varbin}
\alias{print.varbin}
\encoding{latin1}

\title{Estimate of a Probability from Clustered Binomial Data}

\description{The function estimates a probability and its variance from clustered binomial data

\{\eqn{(n_1, m_1), (n_2, m_2), ..., (n_N, m_N)}\},

where \eqn{n_i} is the size of cluster \eqn{i}, \eqn{m_i} the number of \dQuote{successes} (proportions are \eqn{y = m/n}), and \eqn{N} the number of clusters. Confidence intervals are calculated using a normal approximation, which might be inappropriate when the probability is close to 0 or 1.}

\usage{
  varbin(n, m, alpha = 0.05, R = 5000)
  
  \method{print}{varbin}(x, ...)
  }

\arguments{
  \item{n}{A vector of the sizes of the clusters.}
  \item{m}{A vector of the numbers of successes (proportions are \eqn{y = m / n}).}
  \item{alpha}{The significance level for the confidence intervals. Default to 0.05, providing 95\% CI's.}
  \item{R}{The number of bootstrap replicates to compute bootstrap mean and variance. Default to 5000.}
  \item{x}{An object of class \dQuote{varbin}.}
  \item{...}{Further arguments to be passed to \dQuote{print}.}
  }

\details{
  Five methods are used for the estimations. Let us consider \eqn{N} clusters of sizes \eqn{n_1, \ldots, n_N} with observed count responses \eqn{m_1, \ldots, m_N}. We note \eqn{y_i = m_i/n_i (i = 1, \ldots, N)} the observed proportions. The underlying assumption is that the probability, say \eqn{mu}, is homogeneous across the clusters.
  
  \bold{Binomial method:} the probability estimate and its variance are calculated by
  
  \eqn{\mu = (sum_{i} (m_i)) / (sum_{i} (n_i))} (ratio estimate) and
  
  \eqn{\mu * (1 - \mu) / (sum_{i} (n_i) - 1)}, respectively.
  
  \bold{Ratio method:} the probability \eqn{\mu} is estimated as for the binomial method (ratio estimate). The one-stage cluster sampling formula is used to calculate the variance of \eqn{\mu} (see Cochran, 1999, p. 32 and p. 66).
  
  \bold{Arithmetic method:} the probability is estimated by \eqn{\mu = sum_{i} (y_i) / N}. The variance of \eqn{\mu} is estimated by \eqn{sum_{i} (y_i - \mu)^2 / (N * (N - 1))}.

   \bold{Jackknife method:} the probability is estimated by \eqn{\mu} defined by the arithmetic mean of the pseudovalues \eqn{y_{v,i}}. The variance is estimated  by \eqn{sum_{i} (y_{v,i} - \mu)^2 / (N * (N - 1))} (Gladen, 1977, Paul, 1982).
   
    \bold{Bootstrap method:} \eqn{R} samples of clusters of size \eqn{N} are drawn with equal probability from the initial sample \eqn{(y_1, \ldots , y_N)} (Efron and Tibshirani, 1993). The bootstrap estimate \eqn{\mu} and its estimated variance  are the arithmetic mean and the empirical variance (computed with denominator \eqn{R - 1}) of the \eqn{R} binomial ratio estimates, respectively.
  }

\value{An object of class \code{varbin}, printed with \code{print.varbin}.}

\references{
  Cochran, W.G., 1999, 3th ed. \emph{Sampling techniques}. Wiley, New York.\cr
  Efron, B., Tibshirani, R., 1993. \emph{An introduction to the bootstrap}. Chapman and Hall, London.\cr
  Gladen, B., 1977. \emph{The use of the jackknife to estimate proportions from toxicological data in the presence 
  of litter effects}. JASA 74(366), 278-283.\cr
  Paul, S.R., 1982. \emph{Analysis of proportions of affected foetuses in teratological experiments}. 
  Biometrics 38, 361-370.
 }

\seealso{\code{\link[boot:boot]{boot}} }

\examples{
data(rabbits)
z <- rabbits[rabbits$group == "M", ]
varbin(z$n, z$m)
by(rabbits,
	list(group = rabbits$group),
  function(x) varbin(n = x$n, m = x$m, R = 1000))
}

\keyword{datagen}
