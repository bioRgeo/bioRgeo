\name{all_steps}
\alias{all_steps}
\title{Wrap-up function that runs several methods available in the package.}
\usage{
all_steps(dat, site, sp, ab, weight = FALSE)
}
\arguments{
\item{dat}{a data.frame with duplicated rows for sites, one row per species,
a third column indicating the abundance of the species in this site}

\item{input_format}{the format of the input data. It is either a long format data.frame with the occurrences of species per site ('data.frame') or a site-species matrix ('matrix').}

\item{weight}{a logical indicating if abundances should be converted as 0/1.}

\item{site}{the character string corresponding to the site column name.}

\item{sp}{the character string corresponding to the species column name.}

\item{ab}{the character string corresponding to the abundance column name.}

\item{metric}{a character string indicating which similarity metric to chose among \emph{simpson}, \emph{jaccard}, \emph{sorensen}, \emph{whittaker}, \emph{euclidean} and \emph{bray}.}

\item{cluster_method}{Character string to chose among the following options:
         \emph{kmeans}, \emph{meanshift}, \emph{ward.D}, \emph{ward.D2}, \emph{single}, \emph{complete}, \emph{average}, \emph{mcquitty}, \emph{median}, \emph{centroid}, \emph{dbscan}, \emph{gmm}, \emph{diana} or \emph{pam}.}

\item{n_clust}{An integer determining the number of clusters. If set to NULL, an optimal number of clusters will be determined.}

\item{optim_method}{When the chosen method is supervised and there is no a priori regarding the optimal number of clusters, an optimization algorithm is added. The following gap statistic to determine the optimal number of cluster are: \emph{globalmax}, \emph{firstmax}, \emph{Tibs2001SEmax}, \emph{firstSEmax} or \emph{globalSEmax}.}

\item{nstart}{An integer determining the number of random centroids to start K-means analysis and to find the optimal number of clusters.}

\item{B}{An integer determining the number of Monte Carlo bootstrap samples. Necessary when an  algorithm is run to determine the optimal number of clusters.}

\item{K.max}{A numeric determining the maximum number of clusters to consider. Necessary when an  algorithm is run to determine the optimal number of clusters.}
}

\item{network_method}{An algorithm used to perfom the bipartite community detection.
Available options are \emph{greedy}, \emph{betweenness}, \emph{walktrap}, \emph{louvain}, \emph{infomap}, \emph{spinglass}, \emph{leading_eigen}, \emph{label_prop}, \emph{netcarto}, \emph{oslom}, \emph{beckett} and \emph{quanbimo}.}

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

}
\description{
Wrap-up function that runs several methods available in the package.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
sp = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

all_steps(
  dat = ex,
  input_format = "data.frame", weight = FALSE, site = "site", sp = "sp",
  ab = NULL, metric = "simpson",
  cluster_method = c("kmeans", "diana"),
  optim_method = "firstSEmax", n_clust = 2, nstart = 25, B = 50, K.max = 2,
  network_method = c("louvain", "greedy"),
  N = 10, n_runs = 10, t_param = 0.1, cp_param = 0.5, hr = 0)

}
