\name{cluster}
\alias{cluster}
\title{Clustering analysis}

\description{
Assign the sites of a contingency matrix to clusters.
}

\usage{
cluster(dat, method = "ward.D2", n_clust = 4)
}

\arguments{
\item{dat}{a contingency matrix}

\item{method}{Character string to chose among the following options:
         'kmeans', 'meanshift', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid', 'dbscan', 'gmm', 'diana' or 'pam'.}

\item{n_clust}{An integer determining the number of clusters. If set to NULL, an optimal number of clusters will be determined.}

\item{optim_method}{When the chosen method is supervised and there is no a priori regarding the optimal number of clusters, an optimization algorithm is added. The following gap statistic to determine the optimal number of cluster are: 'globalmax', 'firstmax', 'Tibs2001SEmax', 'firstSEmax' or 'globalSEmax'.}

\item{nstart}{An integer determining the number of random centroids to start K-means analysis and to find the optimal number of clusters.}

\item{B}{An integer determining the number of Monte Carlo bootstrap samples. Necessary when an  algorithm is run to determine the optimal number of clusters.}

\item{K.max}{A numeric determining the maximum number of clusters to consider. Necessary when an  algorithm is run to determine the optimal number of clusters..}
}

\details{
The acronym given as argument to the \code{cluster} function correspond to the following methods:
\code{"kmeans"} => K-Means

\code{"pam"} => Partitioning Around Medoids

\code{"meanshift"} => Mean-Shift

\code{"dbscan"} => Density-Based Spatial Clustering of Applications with Noise

\code{"gmm"} => Gaussian Mixture Models

\code{"ward.D"} => Ward.D

\code{"ward.D2"} => Ward.D2

\code{"single"} => single

\code{"complete"} => complete

\code{"average"} => Unweighted pair-group method using arithmetic averages (UPGMA, average linkage)

\code{"mcquitty"} => Weighted pair-group method using arithmetic averages (WPGMA, McQuitty’s method)

\code{"median"} => Weighted pair-group method using centroids (WPGMC, median linkage)

\code{"centroid"} => Unweighted pair-group method using centroids (UPGMC, centroid linkage)

\code{"diana"} = > DIvisive ANAlysis Clustering

These algorithms were implemented in other R packages, which \code{cluster} call.
See \code{\link[stats:hclust]{hclust}} from \pkg{\link{stats}}, \code{\link[cluster:diana]{diana}} from \pkg{\link{cluster}}, \code{\link[fastcluster:hclust]{hclust}} from \pkg{\link{fastcluster}}, \code{\link[mclust:Mclust]{Mclust}} from \pkg{\link{mclust}}, \code{\link[dbscan:dbscan]{dbscan}} from \pkg{\link{dbscan}} and \code{\link[meanShiftR:meanShift]{meanShift}} from \pkg{\link{meanShiftR}} for further details about the available methods.

}

\value{
A data frame with two columns: the sites (or rows of the contingency table) and the clusters assigned to each site for the chosen method.
}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

tmp <- contingency(dat = ex, sp = "species", site = "sites", ab = "count")
cluster(tmp, method = "ward.D2", optim_method = "firstSEmax",
                         nstart = 25, B = 50, K.max = 2)
}
