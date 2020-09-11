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
Available options are \emph{greedy}, \emph{girvan}, \emph{walktrap}, \emph{louvain}, \emph{infomap}, \emph{infomap}, \emph{spinglass}, \emph{leading_eigen}, \emph{label_prop}, \emph{optimal}, \emph{beckett}, \emph{dorman} and \emph{netcarto}.}

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
Run a community detection algorithm on a \code{bipartite} network.
Algorithms available are from igraph package (\emph{greedy}, \emph{girvan}, \emph{walktrap}, \emph{louvain}, \emph{infomap}, \emph{spinglass}, \emph{leading_eigen}, \emph{label_prop} or \emph{optimal}), from the \code{bipartite} package (\emph{beckett} and \emph{dorman}) and from rnetcarto R package (\code{netcarto}).
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
                         algo = "greedy", saving_directory = NULL,
                         weight = FALSE, input = "data frame", site = "site",
                         sp = "sp", ab = NULL)
}
