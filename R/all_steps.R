
all_steps <- function(
  dat, input_format = "data.frame", weight = FALSE,
  site, sp, ab = NULL,
  metric = "simpson",
  cluster_method = NULL,
  optim_method = "firstSEmax", n_clust = NULL, nstart = 25, B = 50, K.max = 20,
  network_method = NULL,
  N = 10, n_runs = 10, t_param = 0.1, cp_param = 0.5, hr = 0){

  ## Controls ----
  # Packages
  require(dplyr)
  require(tidyr)

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
    if(network_method %in% c("beckett", "quanbimo") &
       input_format == "data.frame"){
      warning("The chosen algorithm needs a site-species matrix as input.
            We here perform bioRgeo::contingency() function to convert it but
            you can do this step a priori to save time.")
      dat <-  contingency(dat, site = site, sp = sp, ab = ab, weight = weight)
    }

    if(network_method %in% c("netcarto", "quanbimo")){
      warning("The chosen algorithm is slow. Depending on the size of the
            network, it may take a long time to run.")
    }

    if(!(abs(N - round(N)) < .Machine$double.eps^0.5)){
      stop("N must be an integer setting the number of outer-most loops to run
         before picking the best solution.")
    }

    # Controls for OSLOM parameters
    if(network_method == "oslom"){
      # Controls: dat must be a data.frame with three columns containing id1,
      # id2 and similarity metric
      if(!(is.data.frame(dat))){
        stop("dat must be a data.frame with three columns containing id1, id2
        and a similarity metric.
        It is the output of 'bioRgeo::simil()' function")
      }

      if(ncol(dat) != 3){
        stop("dat must be a data.frame with three columns containing id1, id2
  and similarity metric")
      }

      if(!(is.numeric(dat[, 3]))){
        stop("dat must be a data.frame with three columns containing id1, id2
  and similarity metric")
      }

      if(0 %in% dat[, 3]){
        stop("OSLOM needs strictly positive weights to run. Remove the useless
         lines from the input data.frame.")
      }

      if(!(abs(n_runs - round(n_runs)) < .Machine$double.eps^0.5)){
        stop("n_runs must be an integer setting the number of runs.")
      }

      if(!(is.numeric(t_param))){
        stop("t_param must be numeric.")
      }

      if(t_param > 1 | t_param < 0){
        stop("t_param must be comprised between 0 and 1.")
      }

      if(!(is.numeric(cp_param))){
        stop("cp_param must be numeric.")
      }

      if(cp_param > 1 | cp_param < 0){
        stop("cp_param must be comprised between 0 and 1.")
      }

      if(!(abs(hr - round(hr)) < .Machine$double.eps^0.5)){
        stop("hr must be an integer setting the number of hierarchical levels.")
      }

      if(hr < 0){
        stop("hr must be positive.")
      }
    }

  }

  ## contingency() ----
  # Create contingency table if needed
  if(input_format == "data.frame"){
    sp_mat <- contingency(dat, site = site, sp = sp, ab = ab, weight = weight)
  } else if(input_format == "matrix"){
    sp_df <- contingency_to_df(dat, col1 = "site", col2 = "sp", col3 = "ab",
                               remove_zeros = TRUE)
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
  if(input == "data.frame"){
    bioregions <- dat %>%
      distinct(site, .keep_all = TRUE) # remove duplicates per site
  }else if(input_format == "matrix"){
    bioregions <- sp_df %>%
      distinct(site, .keep_all = TRUE)
  }

  bioregions <- bioregions %>%
    left_join(km_res, by = "site") %>% # add kmeans
    rename(cluster_km = cluster) %>%
    left_join(mshift_res, by = "site") %>% # add meanshift
    rename(cluster_mshift = cluster) %>%
    left_join(oslom_mod, by = "site") %>% # add oslom
    rename(cluster_oslom = bioregion) %>%
    left_join(ward_res, by = "site") %>% # add ward
    rename(cluster_ward = cluster) %>%
    left_join(dbscan_res, by = "site") %>% # add dbscan
    rename(cluster_dbscan = cluster) %>%
    left_join(gmm_res, by = "site") %>% # add gmm
    rename(cluster_gmm = cluster) %>%
    left_join(pam_res, by = "site") %>% # add pam
    rename(cluster_pam = cluster) %>%
    left_join(diana_res, by = "site") %>% # add diana
    rename(cluster_diana = cluster) %>%
    left_join(bip_site, by = "site") %>% # add fastgreedy
    rename(cluster_fastgreedy = module) %>%
    left_join(bip_site_beckett, by = "site") %>% # add lpawb
    rename(cluster_beckett = module) %>%
    left_join(infomap_site, by = "site") %>% # add infomap
    rename(cluster_infomap = module) %>%
    left_join(bip_site_louvain, by = "site") %>% # add louvain
    rename(cluster_louvain = module) %>%
    left_join(bip_site_walktrap, by = "site") %>% # add walktrap
    rename(cluster_walktrap = module) %>%
    as.data.frame()

  all100 <- comparison(bioregions, site = "site",
                       bio_col = c(6:ncol(bioregions)), output = "both")

  ## contribute() ----
  tmp <- left_join(dat, list_res["oslom"], by = "site")
  scores <- contribute(dat = tmp, sp_col = "sp", site_col = "site",
                       bioregion_col = "bioregion")

  ## interact() ----
  lambda <- interact(input_network = "projected",
                     dat = scores, plot = TRUE, output_format = "matrix")

  ## Return outputs ----
  if(!is.null(cluster_method)){
    if(!is.null(network_method)){
      return(list(dat, dat_proj, res_cluster, res_network, all100))
    } else{
      return(list(dat, dat_proj, res_cluster, all100))
    }
  } else{
    return(list(dat, dat_proj, res_network, all100))
  }
}
