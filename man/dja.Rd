\name{dja}
\alias{dja}
\docType{data}
\encoding{latin1}

\title{ Mortality of Djallonke Lambs in Senegal }

\description{
  Field trial to assess the effect of ewes deworming (prevention of gastro-intestinal parasitism) on the mortality
  of their offspring (age < 1 year). This data set is extracted from a large database on small ruminants production
  and health in Senegal (Lancelot et al., 1998). Data were collected in a sample of herds in Kolda (Upper Casamance,
  Senegal) during a multi-site survey (Faug?re et al., 1992). See also the references below for a presentation of the
  follow-up survey (Faug?re and Faug?re, 1986) and a description of the farming systems (Faug?re et al., 1990).
  }

\usage{data(dja)}

\format{
  A data frame with 21 observations on the following 4 variables.
  \describe{
    \item{group}{a factor with 2 levels: \code{CTRL} and \code{TREAT}, indicating the treatment.}
    \item{village}{a factor indicating the village of the herd.}
    \item{herd}{a factor indicating the herd.}
    \item{n}{a numeric vector: the number of animals exposed to mortality.}
    \item{trisk}{a numeric vector: the exposition time to mortality (in year).}
    \item{m}{a numeric vector: the number of deaths.}
  }
}

\references{
  Faug?re, O., Faug?re, B., 1986. \emph{Suivi de troupeaux et contr?le des performances individuelles des petits ruminants en milieu traditionnel africain. Aspects m?thodologiques.} Rev. Elev. M?d. v?t. Pays trop., 39 (1): 29-40.\cr
  Faug?re, O., Dock?s, A.-C., Perrot, C., Faug?re, B., 1990. \emph{L'?levage traditionnel des petits ruminants
  au S?n?gal. I. Pratiques de conduite et d'exploitation des animaux chez les ?leveurs de la r?gion de Kolda.} Revue
  Elev. M?d. v?t. Pays trop. 43: 249-259.\cr
  Faug?re, O., Tillard, E., Faug?re, B., 1992. \emph{Prophylaxie chez les petits ruminants au S?n?gal : r?gionalisation d'une politique nationale de protection sanitaire}. In: B. Rey, S. H. B. Lebbie, L. Reynolds (Ed.), First biennial conference of the African Small Ruminant Research Network, ILCA, 1990, ILRAD, Nairobi, pp. 307-314.\cr
  Lancelot, R., Faye, B., Juan?s, X., Ndiaye, M., P?rochon, L., Tillard, E., 1998. \emph{La base de donn?es BAOBAB: un outil pour mod?liser la production et la sant? des petits ruminants dans les syst?mes d'?levage traditionnels au S?n?gal.} Revue Elev. M?d. v?t. Pays trop., 51 (2): 135-146.\cr
  }

\keyword{datasets}
