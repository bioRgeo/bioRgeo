\name{contingency_to_df}
\alias{contingency_to_df}
\title{Create a data frame from a contingency table.}
\usage{
contingency_to_df(dat, site.col = "site", species.col = "sp", abundance.col = "ab", remove_absent_species = TRUE)
}
\arguments{
\item{dat}{a contingency matrix with sites as rows and species as columns.}

\item{site.col}{character string describing the rows of the contingency matrix.}

\item{species.col}{character string describing the columns of the contingency matrix.}

\item{abundance.col}{character string describing the values of the contingency matrix.}

\item{remove_absent_species}{boolean determining whether absent
    species from the contingency matrix have to be removed from the output.}

}
\description{
Converts a contingency table into a data.frame with replicated rows per site.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

mat <- contingency(dat = ex, site = "sites", sp = "species", ab = "count")
contingency_to_df(mat, site.col = "site", species.col = "sp", abundance.col = "ab")
}
