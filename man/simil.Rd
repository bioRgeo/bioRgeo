\name{simil}
\alias{simil}
\title{Similarity matrix}
\description{
Project a contingency matrix into a similarity data.frame, the similarity being based upon the sites each pair of species has in common.
}
\usage{
simil(dat, metric = "simpson", input = "matrix", output = "data frame", site = NULL, sp = NULL, ab = NULL, weight = FALSE)
}
\arguments{
\item{dat}{a contingency matrix with sites as rows and species as columns.}

\item{metric}{a character string indicating which similarity metric to chose among \emph{simpson}, \emph{jaccard}, \emph{sorensen}, \emph{whittaker}, \emph{euclidean} and \emph{bray}.}

\item{input}{a character string equals to 'matrix' or 'data frame' indicating whether the input of the function is long-format data frame or a wide format matrix.}

\item{output}{a character string indicating whether the output of the function should be a symmetrical matrix, a long-format data frame or a distance object.}

\item{site}{If \code{input = 'data frame'}, the column name of the sites.}

\item{sp}{If \code{input = 'data frame'}, the column name of the species.}

\item{ab}{If \code{input = 'data frame'} and species have abundances, the column name of the abundances of species per site.}

\item{weight}{a boolean indicating whether abundances should be taken into account.}
}
\details{
With \code{a} the number of species shared by a pair of sites, \code{b} species only present in the first site and \code{c} species only present in the second site.

\eqn{Simpson = 1 - min(b, c)/(a + min(b, c))}

\eqn{Sorensen = 1 - (b+c)/(2a + b + c)}

\eqn{Jaccard = 1 - a/(a + b + c)}

\eqn{Whittaker = 1 - (a + b + c)/((2a + b + c)/2)}

Euclidean computes the Euclidean distance between each pair of site.

If abundances data are available, Bray-Curtis dissimilarity can also be computed. It computes following this equation:
\eqn{Bray-Curtis = 1 - 2C_{ij}/(S_i + S_j)}
with \eqn{2C_{ij}} the sum of the lesser values for common species shared by a pair of sites. \eqn{S_{i}} and \eqn{S_{j}} are the total number of specimens counted at both sites \eqn{i} and \eqn{j}.

}
\value{
If \code{output = 'data frame'}, a data.frame with 3 columns. The first two columns, \code{id1} and \code{id2}, correspond to each pair of sites. The third column indicates the similarity value for the chosen metric between each pair of site.

If \code{output = 'matrix'}, a matrix with n_sites rows and n_sites columns. Values of the matrix indicate the similarity value for the chosen metric between each pair of site.

If \code{output = 'dist'}, a dist object based on the symmetrical matrix.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

ex_mat <- contingency(dat = ex, sp = "species", site = "sites", ab = "count",
weight = FALSE)
simil(ex_mat, metric = "simpson", input = "matrix")
}
