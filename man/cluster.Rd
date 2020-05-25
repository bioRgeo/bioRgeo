\name{cluster}
\alias{cluster}
\title{Clustering analysis from a contingency table}
\usage{
cluster(dat, method = "ward.D2", n_clust = 4)
}
\arguments{
\item{dat}{a contingency table}

\item{method}{Hierarchical clustering method chosen among the followings:
         ward.D, ward.D2, single, complete, average, mcquitty, median,
         centroid, dbscan or gmm.}

\item{optim_method}{Chosen gap statistic to determine the optimal number of
cluster among the followings: globalmax, firstmax, Tibs2001SEmax, firstSEmax
or globalSEmax.}

\item{n_clust}{An integer determining the number of clusters. If set to NULL,
an optimal number of clusters will be determined.}

\item{nstart}{An integer determining the number of random centroids
         to start k-means analysis.}

\item{B}{An integer determining the number of Monte Carlo bootstrap
         samples.}

\item{K.max}{A numeric determining the maximum number of clusters
         to consider.}
}
\description{
Performs a clustering analysis such as k-means clustering, determines the
optimal number of clusters and extract them.
}
\details{}
\section{Note}{}

\arguments{
\item{}{a long format data.frame with sites as rows and the associated
cluster.}
}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

tmp <- contingency(dat = ex, sp = "species", site = "sites", ab = "count")
cluster(tmp, method = "ward.D2", optim_method = "firstSEmax",
                         nstart = 25, B = 50, K.max = 2)
}
