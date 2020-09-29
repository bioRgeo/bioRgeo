
all_steps <- function(
  dat, input_format = "data.frame", weight = FALSE,
  site, sp, ab = NULL,
  metric = "simpson",
  cluster_method = NULL,
  optim_method = "firstSEmax", n_clust = NULL, nstart = 25, B = 50, K.max = 20,
  network_method = NULL,
  weight = FALSE, N = 10, n_runs = 10, t_param = 0.1, cp_param = 0.5, hr = 0){

  ## Controls ----
  if(!is.character(input_format)){
    stop("input_format designs the format of the input data. It is either a
    long format data.frame with the occurrences of species per site
    ('data.frame') or a site-species matrix ('matrix').")
  }

  if(!(input_format %in% c("data.frame", "matrix"))){
    stop("input_format designs the format of the input data. It is either a
    long format data.frame with the occurrences of species per site
    ('data.frame') or a site-species matrix ('matrix').")
  }

  if(input_format == "matrix"){
    if(!is.matrix(dat)){
      stop("dat should be a matrix with sites in rows and species in columns.")
    }

    if(is.null(rownames(dat)) | is.null(colnames(dat))){
      stop("dat does not have rownames or colnames.")
    }
  }

  if(input_format == "data.frame"){
    if(!is.data.frame(dat)){
      stop("dat must be a long format data.frame with the occurrences of
           species per site. It must have a column for sites and another
           one for species. A third optional column for abundances of
           species within sites can be present.")
    }

    if(!is.character(site)){
      stop("site must be the column name of dat describing the sites.")
    }

    if(!is.character(sp)){
      stop("sp must be the column name of dat describing the occurrences
           of species within sites.")
    }

    if(!is.null(ab) & !is.character(ab)){
      stop("ab must be the column name of dat describing the abundances
         of species within sites.")
    }

    if(is.null(ab) & weight == TRUE){
      warning("Without column abundances, contingency table will only get binary
        values.")
    }
  }

  if(!(metric %in% c("simpson", "jaccard", "sorensen", "whittaker", "bray",
                     "euclidean"))){
    stop("Similarity metric chosen is not available.
     Please chose among the followings:
         'simpson', 'jaccard', 'sorensen', 'whittaker', 'bray' or
         'euclidean.'")
  }

  if(is.null(cluster_method) & is.null(network_method)){
    stop("You should provide at least one bioregionalization method, either
         clustering or network algorithm.")

  }

  if(!is.null(cluster_method)){
    if(!(cluster_method %in% c("kmeans", "meanshift",
                               "ward.D", "ward.D2", "single", "complete",
                               "average", "mcquitty", "median", "centroid", "dbscan",
                               "gmm", "diana", "pam"))){
      stop("Clustering algorithm chosen is not available.
     Please chose among the following clustering techniques:
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
    }else if(class(dat) == "data.frame"){
      if(K.max > length(unique(dat[, site]))){
        stop("K.max should not be superior to the number of sites.")
      }
    }
  }

  if(!is.null(network_method)){
    if(!(network_method %in% c("greedy", "betweenness", "walktrap", "louvain",
                               "LPAwb", "infomap", "spinglass", "leading_eigen",
                               "label_prop", "netcarto", "oslom", "infomap",
                               "beckett", "quanbimo"))){
      stop("Network algorithm chosen is not available.
     Please chose among the following network algorithms:
         'greedy', 'betweenness', 'walktrap', 'louvain', 'LPAwb',
         'spinglass', 'leading_eigen', 'label_prop', 'netcarto', 'oslom',
         'infomap', 'beckett' or 'quanbimo.'")
    }
  }

  ## contingency() ----
  # Create contingency table if needed
  if(input_format == "data.frame"){
    sp_mat <- contingency(dat, site = site, sp = sp, ab = ab, weight = weight)
  }

  ## simil() ----
  # Project network with simil()
  dat_proj <- simil(sp_mat, metric = metric, input = "matrix",
                    output = "data frame",
                    site = NULL, sp = NULL, ab = NULL, weight = FALSE)

  dat_proj <- dat_proj[, c("id1", "id2", metric)]

  ## cluster() ----
  if(!is.null(cluster_method)){
    if(input_format == "matrix"){
      res_cluster <- cluster(dat, method = cluster_method,
                             optim_method = optim_method, n_clust = n_clust,
                             nstart = nstart, B = B, K.max = K.max)
    } else if(input_format == "data.frame"){
      res_cluster <- cluster(sp_mat, method = cluster_method,
                             optim_method = optim_method, n_clust = n_clust,
                             nstart = nstart, B = B, K.max = K.max)
    }
  }

  ## community() ----
  if(!is.null(network_method)){
    res_network <- community(dat = sp_mat, algo = network_method,
                             input = input_format,
                             N = N, n_runs = n_runs, t_param = t_param,
                             cp_param = cp_param, hr = hr, weight = weight)
  }

  ## comparison() ----
  # Gather all the bioregionalizations
  all_bioregions <- dat %>%
    select(site, x, y, env) %>%
    distinct(site, .keep_all = TRUE) %>%
    left_join(list_res["oslom"], by = "site") %>% # add OSLOM
    rename(oslom = bioregion) %>%
    left_join(list_res["ward"], by = "site") %>% # add Ward
    rename(ward = cluster)

  # Test of comparison function
  all100 <- comparison(all_bioregions, bio_col = c(5:6))

  ## contribute() ----
  tmp <- left_join(dat, list_res["oslom"], by = "site")
  scores <- contribute(dat = tmp, sp_col = "sp", site_col = "site",
                       bioregion_col = "bioregion")

  ## interact() ----
  lambda <- interact(input_network = "projected",
                     dat = scores, plot = TRUE, output_format = "matrix")

  ## Return outputs ----
  return(list(dat, dat_proj))

}
