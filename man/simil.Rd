\name{simil}
\alias{simil}
\title{Similarity matrix}
\description{
Project a contingency matrix into a similarity data.frame, the similarity being based upon the sites each pair of species has in common.
}
\usage{
simil(contingency_mat, metric = "simpson")
}
\arguments{
\item{contingency_mat}{a contingency matrix with sites as rows and species as columns.}

\item{metric}{a character string indicating which similarity metric to chose among "simpson", "jaccard", "sorensen", "whittaker" and "euclidean".}

}
\details{
With \code{a} the number of species shared by a pair of sites, \code{b} species only present in the first site and \code{c} species only present in the second site.

\eqn{Simpson = 1 - min(b, c)/(a + min(b, c))}

\eqn{Sorensen = 1 - (b+c)/(2a + b + c)}

\eqn{Jaccard = 1 - a/(a + b + c)}

\eqn{Whittaker = 1 - (a + b + c)/((2a + b + c)/2)}

Euclidean computes the Euclidean distance between each pair of site.
}
\value{
A data.frame with 3 columns. The first two columns, \code{id1} and \code{id2}, correspond to each pair of sites. The third column indicates the similarity value for the chosen metric between each pair of site.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

ex_mat <- contingency(dat = ex, sp = "species", site = "sites", ab = "count",
binary = TRUE)
simil(ex_mat, metric = "simpson")
}
