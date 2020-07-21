
cluster <- function(dat, method = "ward.D2", optim_method = "firstSEmax",
                    n_clust = NULL, nstart = 25, B = 50, K.max = 20){
  require(cluster)

  if(!is.data.frame(dat)){
    stop("dat should be the output of the similarity function, i.e. a data
         frame with three columns. The first two columns contains each pair of
         sites and the third column their similarity value.")
  }

  if(!(method %in% c("ward.D", "ward.D2", "single", "complete", "average",
                     "mcquitty", "median", "centroid", "dbscan", "gmm"))){
    stop("Hierarchical clustering method chosen is not available.
     Please chose among the followings:
         ward.D, ward.D2, single, complete, average,
         mcquitty, median, centroid, dbscan or gmm.")
  }

  if(!(optim_method %in% c("globalmax", "firstmax", "Tibs2001SEmax",
                           "firstSEmax", "globalSEmax."))){
    stop("Chosen gap statistic to determine the optimal number of cluster is
    not available.
     Please chose among the followings:
         globalmax, firstmax, Tibs2001SEmax, firstSEmax or globalSEmax.")
  }

  if(!is.null(n_clust)){
    if(!(abs(n_clust - round(n_clust)) < .Machine$double.eps^0.5)){
      stop("n_clust must be an integer determining the number of clusters.")
    }
  }

  if(!(abs(nstart - round(nstart)) < .Machine$double.eps^0.5)){
    stop("nstart must be an integer determining the number of random centroids
         to start k-means analysis.")
  }

  if(!(abs(B - round(B)) < .Machine$double.eps^0.5)){
    stop("B must be an integer determining the number of Monte Carlo bootstrap
         samples.")
  }

  if(!is.numeric(K.max)){
    stop("K.max must be a numeric determining the maximum number of clusters
         to consider.")
  }

  if(K.max > length(unique(dat[, 1]))){
    stop("K.max should not be superior to the number of sites.")
  }

  if(method == "dbscan"){
    require(dbscan)
    # Size of the epsilon neighborhood: normally determined by observing the
    # knee-plot
    eps <- quantile(kNNdist(dat, k = 5)[
      kNNdist(dat, k = 5) > 0], 0.25)
    # eps <- mean(kNNdist(dat, k = 5))
    db_clust <- dbscan(x = dat, eps = eps, minPts = 5)

    # Data.frame of results
    res <- data.frame(site = rownames(dat),
                      cluster = as.character(db_clust$cluster))

  } else if(method == "gmm"){
    # Conversion to dist object
    dat <- as.dist(dat)
    require(mclust)
    gmm_mclust <- Mclust(dat)

    # Data.frame of results
    res <- data.frame(site = names(gmm_mclust$classification),
                      cluster = as.character(gmm_mclust$classification))

  } else{
    # Conversion to dist object
    dat <- as.dist(dat)

    #h <- hclust(dat, method = method)
    require(fastcluster)
    h <- fastcluster::hclust(dat, method = method)
    # plot(h)

    if(!is.null(n_clust)){
      optim_k <- n_clust
    } else{

      # Determine optimal numbers of clusters
      # https://stackoverflow.com/questions/53159033/how-to-get-the-optimal-number-of-clusters-from-the-clusgap-function-as-an-output
      # https://uc-r.github.io/kmeans_clustering#gap
      # https://uc-r.github.io/hc_clustering

      gap_stat <- clusGap(dat, FUN = kmeans, nstart = nstart, K.max = K.max,
                          B = B)

      optim_k <- maxSE(f = gap_stat$Tab[, "gap"],
                       SE.f = gap_stat$Tab[, "SE.sim"],
                       method = optim_method)
    }

    # Cut the tree with optim_k numbers
    dend <- as.dendrogram(h)
    # Data.frame of results
    res <- data.frame(site = names(dendextend::cutree(dend,
                                                      k = optim_k)),
                      cluster = as.character(dendextend::cutree(dend,
                                                                k = optim_k)))
  }
  return(res)
  # Visualization
  # factoextra::fviz_nbclust(dat, kmeans, method = "gap_stat", k.max = 20)

  # dend <- as.dendrogram(h)
  # library(dendextend)
  # labels(dend) <- x[order.dendrogram(dend)]
  # # due to the ties - there is specific reason to have this be these 3 clusters:
  # cutree(dend, k = 3)
  #
  # library(NbClust)
  # NbClust(data = dat, diss = NULL, distance = "euclidean",
  #         min.nc = 1, max.nc = 15, method = "ward.D2")
}

