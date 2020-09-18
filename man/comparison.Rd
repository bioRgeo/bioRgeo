\name{comparison}
\alias{comparison}
\title{Function to identify pixels grouped together through different
bioregionalizations.}
\usage{
comparison(dat, bio_col, thres = 10, output = "both")
}
\arguments{
\item{dat}{a data.frame with a column containing the sites and
columns containing different partitions.}

\item{bio_col}{a numeric indicating which columns of dat contains
the different bioregions.}

\item{thres}{a numeric determining the minimum number of pixels per pair.}

\item{output}{A character string indicating whether the output format is a square matrix with the number of times each pair of sites have been classed together, a data frame with the sites all the time classed together or both objects.}

}
\description{
Function to identify pixels grouped together through different methods of
bioregionalizations.

}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(site = paste0("site", seq(1:7)),
bio1 = c("a", "b", "a", "c", "d", "b", "d"),
bio2 = c("aa", "aa", "aa", "bb", "cc", "bb", "cc"))

comparison(ex, site = "site", bio_col = c(2, 3), thres = 0, output = "both")
}
