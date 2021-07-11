\name{interact}
\alias{interact}
\title{Computes interactions between modules.}
\usage{
interact(dat, site, sp, bio_site, plot = FALSE, chord = FALSE)
}
\arguments{
\item{dat}{a data.frame with four columns containing the sites, the species occurrences in each site and the bioregion of each site.}

\item{site}{a character string corresponding to the column with sites.}

\item{sp}{a character string corresponding to the column with species occurrences within sites.}

\item{bio_site}{a character string corresponding to the column with the bioregion of each site.}

\item{plot}{a boolean indicating whether to store the heatmap plot.}

\item{chord}{a boolean indicating whether to plot the chord diagram.}
}
\description{
See Lenormand et al., Ecology and Evolution, 2019.
}
\details{
compute a Sorensen index for each pair of bioregions according the following formula:
\deqn{
\beta_{ij, i \neq j} = \frac{2a}{b + c}
}

With \eqn{i} and \eqn{j} two distinct bioregions, \eqn{a} the number of species commonly occurring in their sites, \eqn{b} the number of species occurring in bioregion \eqn{i} and \eqn{c} the number of species occurring in bioregion \eqn{j}. This index ranks between 0, when no species are shared between two bioregions, and 1, when the species list are identical.
A heatmap and a chord diagram can be plotted.
}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20),
bioregion = c(rep("region1", 5), rep("region2", 2)))

interact(dat = ex, site = "sites", sp = "species", bio_site = "bioregion", plot = TRUE, chord = TRUE)
}
