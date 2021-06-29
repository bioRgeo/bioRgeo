\name{contingency}
\alias{contingency}
\title{Create a contingency table from a data.frame}
\usage{
contingency(dat, site.col, species.col, abundance.col)
}
\arguments{
\item{dat}{a data.frame with duplicated rows for sites, one row per species,
a third column indicating the abundance of the species in this site}

\item{site.col}{the character string corresponding to the site column name}

\item{species.col}{the character string corresponding to the species column name}

\item{abundance.col}{the character string corresponding to the abundance column name}

}
\description{
Converts a data.frame into a contingency table.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

contingency(dat = ex, site.col = "sites", species.col = "species", abundance.col = "count")
}
