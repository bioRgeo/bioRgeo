\name{community}
\alias{community}
\title{Community detection on bipartite networks}
\usage{
community(dat, algorithm = "greedy", weight = FALSE, input = "matrix", site.col = NULL, species.col = NULL, abundance.col = NULL)
}
\arguments{
\item{dat}{a contingency matrix with sites as rows and species as columns.If 'algorithm="oslom"', a data.frame with three columns containing id1, id2 and a similarity metric (output of 'bioRgeo::simil()' function).}

\item{algorithm}{An algorithm used to perfom the bipartite community detection.
Available options are \emph{greedy}, \emph{betweenness}, \emph{walktrap}, \emph{louvain}, \emph{infomap}, \emph{spinglass}, \emph{leading_eigen}, \emph{label_prop}, \emph{netcarto}, \emph{oslom}, \emph{beckett} and \emph{quanbimo}.}

\item{weight}{A boolean indicating whether the detection of communities has to
be weighted or not.}

\item{input}{A character string indicating whether the input format is a
contingency matrix or a long-format data.frame.}

\item{site.col}{a character string corresponding to the column with sites.}

\item{species.col}{a character string corresponding to the column with species.}

\item{abundance.col}{a character string corresponding to the column with abundances of
species at sites.}

\item{N}{Integer setting the number of outer-most loops to run before picking the best solution, for Infomap algorithm only.}

\item{n_runs}{Numeric parameter setting the number of runs for the first
hierarchical level. The default value is 10. The bigger, the slower the
function but the more accurate the results.}

\item{t_param}{Numeric parameter setting the probability value to consider
a module as significant. Must be comprised between 0 and 1.}

\item{cp_param}{Numeric parameter setting the resolution. Must be comprised
between 0 and 1.}

\item{hr}{Numeric parameter setting the number of hierarchical levels.
Must be a positive integer.}

\item{oslom_id1}{With oslom algorithm only. A character string corresponding to the column name of the first pair of sites.}

\item{oslom_id2}{With oslom algorithm only. A character string corresponding to the column name of the second pair of sites..}

\item{oslom_proj}{With oslom algorithm only. A character string corresponding to the column name of the similarity metric between each pair of sites.}

}
\description{
Run a community detection algorithm on a bipartite network.\cr
Algorithms available are \emph{greedy}, \emph{betweenness}, \emph{walktrap}, \emph{louvain}, \emph{spinglass}, \emph{leading_eigen}, \emph{label_prop} (all from from \code{igraph} package), \emph{beckett} and \emph{quanbimo} from \code{bipartite} package, \emph{netcarto} from \code{rnetcarto} R package and \emph{oslom} and \emph{infomap} that come from compiled codes (respectively from \code{oslom.org} and \code{mapequation.org/infomap})\cr
Algorithm can be weighted in case of a weighted bipartite network.
}

\value{
A data.frame with 2 columns: the sites and their associated bioregions.\cr
Some algorithms (the ones from \code{igraph} package and \emph{infomap} return a bioregion for both types of nodes from the bipartite network. This means that both sites and species are classified into bioregions.)
}

\details{}
\section{Note}{}

\examples{
# Example dataset from bioRgeo:
data("virtual")

# fastgreedy algorithm from igraph
res <- community(dat = virtual[which(virtual$pa > 0), ],
                         algorithm = "greedy",
                         weight = FALSE, input = "data.frame", site.col = "site",
                         species.col = "sp", abundance.col = NULL)
}
