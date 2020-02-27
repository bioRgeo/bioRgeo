\name{algo_bipartite}
\alias{algo_bipartite}
\title{Community detection on bipartite networks}
\usage{
algo_bipartite(dat, algo = "greedy", weight = FALSE, input = "matrix",
site = NULL, sp = NULL, ab = NULL)
}
\arguments{
\item{dat}{a contingency matrix with sites as rows and species as columns.}

\item{algo}{An algorithm used to perfom the bipartite community detection.
Available options are greedy, girvan, walktrap, louvain, infomap,
infomap, spinglass, leading_eigen, label_prop, optimal, beckett, dorman and
netcarto.}

\item{weight}{A boolean indicating whether the detection of communities has to
be weighted or not.}

\item{input}{A character string indicating whether the input format is a
contingency matrix or a long-format data.frame.}

\item{site}{a character string corresponding to the column with sites.}

\item{sp}{a character string corresponding to the column with species.}

\item{ab}{a character string corresponding to the column with abundances of
species at sites.}

}
\description{
Run a community detection algorithm on a bipartite network.
Algorithms available are from igraph package (greedy, girvan, walktrap,
louvain, infomap, spinglass, leading_eigen, label_prop or optimal),
from the bipartite package (beckett and dorman) and from rnetcarto R package
(netcarto).
Algorithm can be weighted in case of a weighted bipartite network.
}
\details{}
\section{Note}{}

\examples{
# Example dataset from bioRgeo:
data("virtual_sp")
sp_df <- virtual[[1]]

# fastgreedy algorithm from igraph
bip_df <- algo_bipartite(dat = sp_df[which(sp_df$pa > 0), ],
                         algo = "greedy",
                         weight = FALSE, input = "data frame", site = "site",
                         sp = "sp", ab = NULL)
}
