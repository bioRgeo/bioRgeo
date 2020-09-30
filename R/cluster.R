
cluster <- function(dat, method = "ward.D2", optim_method = "firstSEmax",
                    n_clust = NULL, nstart = 25, B = 50, K.max = 20){
  require(cluster)

  ## 1. Controls ----
  if(!is.matrix(dat) & !(class(dat) == "dist")){
    stop("dat should be either a matrix with sites as rows and species as
         columns or a 'dist' object.")
  }

  if(!(method %in% c("kmeans", "meanshift",
                     "ward.D", "ward.D2", "single", "complete",
                     "average", "mcquitty", "median", "centroid", "dbscan",
                     "gmm", "diana", "pam"))){
    stop("Hierarchical clustering method chosen is not available.
     Please chose among the followings:
         'kmeans', 'meanshift',
         'ward.D', 'ward.D2', 'single', 'complete', 'average',
         'mcquitty', 'median', 'centroid', 'dbscan', 'gmm', 'diana' or
         'pam'.")
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

  if(is.null(n_clust) & method %in%
     c("kmeans", "ward.D", "ward.D2", "single", "complete",
       "average", "mcquitty", "median", "centroid", "diana", "pam")){
    warning("The chosen method is a supervised algorithm that needs a number of
            clusters. Since 'n_clust = NULL', an optimization algorithm is
            first executed to determine the optimal numbers of clusters.
            This step may take a while.")
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

  if(is.matrix(dat)){
    if(K.max > nrow(dat)){
      stop("K.max should not be superior to the number of sites.")
    }
  }else if(class(dat) == "dist"){
    if(K.max > length(labels(dat))){
      stop("K.max should not be superior to the number of sites.")
    }
  }
  ## 2. Input conversion ----
  if(!(method %in% c("meanshift", "dbscan", "diana"))){
    if(is.null(n_clust)){
      # Storing the matrix input, necessary to find the optimal nb of clusters
      dat1 <- dat
    }
    if(class(dat) != "dist"){
      # Project dat using simil() function with simpson metric
      dat <- simil(dat, metric = "simpson", output = "dist")
    }
  }

  ## 3. Clustering ----
  if(method == "meanshift"){
    require(meanShiftR)
    tmp <- meanShift(queryData = dat, algorithm = "LINEAR", alpha = 0,
                     iterations = 100)

    # Data.frame of results
    res <- data.frame(site = rownames(dat),
                      cluster = tmp$assignment)

  } else if(method == "dbscan"){
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
    # dat <- as.dist(dat)

    # Number of clusters for supervised algorithms
    if(!is.null(n_clust)){
      optim_k <- n_clust
    } else{
      # Determine optimal numbers of clusters
      # https://stackoverflow.com/questions/53159033/how-to-get-the-optimal-number-of-clusters-from-the-clusgap-function-as-an-output
      # https://uc-r.github.io/kmeans_clustering#gap
      # https://uc-r.github.io/hc_clustering

      if(method == "diana"){
        gap_stat <- clusGap(dat, FUN = kmeans, nstart = nstart, K.max = K.max,
                            B = B)
      } else{
        gap_stat <- clusGap(dat1, FUN = kmeans, nstart = nstart, K.max = K.max,
                            B = B)
      }

      optim_k <- maxSE(f = gap_stat$Tab[, "gap"],
                       SE.f = gap_stat$Tab[, "SE.sim"],
                       method = optim_method)
    }

    if(!(method %in% c("kmeans", "pam"))){
      if(method == "diana"){
        h <- cluster::diana(dat)
      } else{
        # h <- hclust(dat, method = method)
        require(fastcluster)
        h <- fastcluster::hclust(dat, method = method)
        # plot(h)
      }

      # Cut the tree with optim_k numbers
      dend <- as.dendrogram(h)
      # Data.frame of results
      res <- data.frame(site = names(dendextend::cutree(dend,
                                                        k = optim_k)),
                        cluster = as.character(dendextend::cutree(dend,
                                                                  k = optim_k)))
    } else if(method == "kmeans"){
      h <- kmeans(dat, centers = optim_k, iter.max = 10, nstart = nstart)
      res <- data.frame(site = names(h$cluster),
                        cluster = as.numeric(h$cluster))
    } else if(method == "pam"){
      h <- cluster::pam(dat, k = optim_k, metric = "euclidean")
      res <- data.frame(site = names(h$clustering),
                        cluster = as.numeric(h$clustering))
    }
  }

  # Convert cluster column into character
  res$cluster <- as.character(res$cluster)

  # Changing column names of res: paste(cluster, method)
  colnames(res)[colnames(res) == "cluster"] <- paste0("cluster_", method)

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

