
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

  if(!is.logical(weight)){
    stop("'weight' must be a boolean.")
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
    if(!(all(cluster_method %in% c("kmeans", "meanshift",
                                   "ward.D", "ward.D2", "single", "complete",
                                   "average", "mcquitty", "median", "centroid",
                                   "dbscan", "gmm", "diana", "pam")))){
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

    if(is.null(n_clust) & any(
      cluster_method %in%
      c("kmeans", "ward.D", "ward.D2", "single", "complete",
        "average", "mcquitty", "median", "centroid", "diana", "pam"))){
      warning("One of the chosen method is a supervised algorithm that needs a
    number of clusters. Since 'n_clust = NULL', an optimization algorithm is
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

    # if(length(cluster_method) > 1 &
    #    any(c("meanshift", "dbscan", "diana") %in% cluster_method) &
    #    !(all(cluster_method %in% c("meanshift", "dbscan", "diana")))){
    #   stop("'meanshift', 'dbscan' and 'diana' need a different input than the
    # other clustering techniques. You should run cluster() separately with
    #      these methods.")
    # }
  }

  if(!is.null(network_method)){
    if(!(all(network_method %in% c(
      "greedy", "betweenness", "walktrap", "louvain", "infomap", "spinglass",
      "leading_eigen", "label_prop", "netcarto", "oslom", "infomap", "beckett",
      "quanbimo")))){
      stop("Provided algorithm to detect communities is not available. Please
    choose among the followings: 'greedy', 'betweenness', 'walktrap', 'louvain',
         'spinglass', 'leading_eigen', 'label_prop', 'netcarto', 'oslom',
         'infomap', 'beckett' or 'quanbimo'")
    }

    if(any(network_method %in% c("beckett", "quanbimo")) &
       input_format == "data.frame"){
      warning("beckett and quanbimo algorithms need a site-species matrix as
    an input. We here perform bioRgeo::contingency() function to convert it but
        you can do this step a priori to save time.")
      dat_to_matrix <- contingency(dat, site = site, sp = sp, ab = ab,
                                   weight = weight)
    }

    if(any(network_method %in% c("netcarto", "quanbimo"))){
      warning("netcarto and quanbimo algorithm is slow. Depending on the size
      of the network, it may take a long time to run.")
    }

    if(!(abs(N - round(N)) < .Machine$double.eps^0.5)){
      stop("N must be an integer setting the number of outer-most loops to run
         before picking the best solution.")
    }

    # Controls for OSLOM parameters
    if("oslom" %in% network_method){
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
        stop("hr must be an integer setting the number of hierarchical
             levels.")
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
                    output = "data.frame",
                    site = NULL, sp = NULL, ab = NULL, weight = FALSE)

  dat_proj <- dat_proj[, c("id1", "id2", metric)]

  ## cluster() ----
  if(!is.null(cluster_method)){
    if(input_format == "matrix"){
      if(length(cluster_method) > 1 &
         any(c("meanshift", "dbscan", "diana") %in% cluster_method) &
         !(all(cluster_method %in% c("meanshift", "dbscan", "diana")))){
        cluster_method1 <- cluster_method[!(cluster_method %in% c(
          "meanshift", "dbscan", "diana"))]
        cluster_method2 <- cluster_method[cluster_method %in% c(
          "meanshift", "dbscan", "diana")]

        res_cluster <- cluster(dat, method = cluster_method1,
                               optim_method = optim_method, n_clust = n_clust,
                               nstart = nstart, B = B, K.max = K.max)

        res_cluster2 <- cluster(dat, method = cluster_method2,
                                optim_method = optim_method, n_clust = n_clust,
                                nstart = nstart, B = B, K.max = K.max)

        # Merge two clustering objects
        res_cluster <- left_join(res_cluster, res_cluster2, by = "site")
        rm(res_cluster2)
      }

    } else if(input_format == "data.frame"){
      if(length(cluster_method) > 1 &
         any(c("meanshift", "dbscan", "diana") %in% cluster_method) &
         !(all(cluster_method %in% c("meanshift", "dbscan", "diana")))){
        cluster_method1 <- cluster_method[!(cluster_method %in% c(
          "meanshift", "dbscan", "diana"))]
        cluster_method2 <- cluster_method[cluster_method %in% c(
          "meanshift", "dbscan", "diana")]

        res_cluster <- cluster(sp_mat, method = cluster_method1,
                               optim_method = optim_method, n_clust = n_clust,
                               nstart = nstart, B = B, K.max = K.max)

        res_cluster2 <- cluster(sp_mat, method = cluster_method2,
                                optim_method = optim_method, n_clust = n_clust,
                                nstart = nstart, B = B, K.max = K.max)

        # Merge two clustering objects
        res_cluster <- left_join(res_cluster, res_cluster2, by = "site")
        rm(res_cluster2)
      }
    }
  }

  ## community() ----
  if(!is.null(network_method)){
    if(length(network_method) > 1 & "oslom" %in% network_method){
      network_method1 <- network_method[!(network_method %in% "oslom")]

      res_network <- community(dat = sp_mat, algo = network_method1,
                               input = "matrix",
                               N = N, n_runs = n_runs, t_param = t_param,
                               cp_param = cp_param, hr = hr, weight = weight)

      res_oslom <- community(dat = dat_proj[which(dat_proj[, metric] > 0), ],
                             algo = "oslom",
                             input = "matrix",
                             N = N, n_runs = n_runs, t_param = t_param,
                             cp_param = cp_param, hr = hr, weight = weight)
    } else{
      res_network <- community(dat = sp_mat, algo = network_method,
                               input = "matrix",
                               N = N, n_runs = n_runs, t_param = t_param,
                               cp_param = cp_param, hr = hr, weight = weight)
    }
  }

  ## comparison() ----
  # Only sites for some network methods
  if(!is.null(network_method)){
    res_network_site <- res_network[which(res_network$cat == "site"), ]
    colnames(res_network_site)[colnames(res_network_site) == "node"] <- "site"
    net_col <- grep("module_", colnames(res_network_site))
    res_network_site <-
      res_network_site[, c("site", colnames(res_network_site)[net_col])]
  }
  # Add oslom if present
  if(!is.null(network_method) & "oslom" %in% network_method){
    res_network_site <- left_join(res_network_site, res_oslom, by = "site")
  }

  # Group all the bioregionalizations in one data.frame
  if(!is.null(cluster_method)){
    if(!is.null(network_method)){
      bioregions <- left_join(res_cluster, res_network_site, by = "site") %>%
        as.data.frame()
    }
    else if(is.null(network_method)){
      bioregions <- res_cluster
    }
  } else if(is.null(cluster_method)){
    bioregions <- res_network_site
  }
  # Delete res_cluster and res_network to save space
  # rm(res_cluster); rm(res_network); rm(res_network_site)

  comp <- comparison(bioregions, site = site,
                     bio_col = c(grep("cluster_", colnames(bioregions)),
                                 grep("module_", colnames(bioregions))),
                     output = "both", thres = 10)

  ## contribute() ----
  # contribute() is run one time per selected method of bioregionalization
  contrib_list <- list()
  interact_list <- list()

  all_method <- c(cluster_method, network_method)
  for(i in 1:length(all_method)){
    # Bioregions for sites only: bioregions

    # Bioregions for species only => some algorithms assign bioregions
    # to species
    if(all_method[i] %in% c("greedy", "betweenness", "walktrap", "louvain",
                            "infomap", "spinglass", "leading_eigen",
                            "label_prop", "netcarto", "infomap", "beckett",
                            "quanbimo")){
      # Merging bioregions of sites with dat
      bioregions_i <- bioregions[, c(1, grep(all_method[i],
                                             colnames(bioregions)))]
      colnames(bioregions_i) <- c("site", "bio_site")
      # Merge with site-species data frame
      if(input_format == "data.frame"){
        dat_i <- left_join(dat, bioregions_i, by = "site")
      } else{
        dat_i <- left_join(sp_df, bioregions_i, by = "site")
      }

      sp_bio_network <- res_network %>%
        filter(cat == "sp") %>%
        rename(sp = node)

      colnames(sp_bio_network)[grep(all_method[i], colnames(sp_bio_network))] <-
        "bio_sp"
      sp_bio_network <- sp_bio_network[, c("sp", "bio_sp")]

      dat_i <- left_join(dat_i, sp_bio_network, by = "sp")

      # assigning species to the bioregion where they occur most for other
      # algorithms
    } else{
      if(input_format == "data.frame"){
        # Merging bioregions of sites with dat
        bioregions_i <- bioregions[, c(1, grep(all_method[i],
                                               colnames(bioregions)))]
        colnames(bioregions_i) <- c("site", "bio_site")
        dat_i <- left_join(dat, bioregions_i, by = "site")

        # Assigning bioregions to species
        sp_bio_network <- dat_i %>%
          group_by(sp, bio_site) %>%
          summarise(c = n()) %>%
          filter(row_number(desc(c)) == 1) %>%
          rename(bio_sp = bio_site) %>%
          select(sp, bio_sp)

        dat_i <- left_join(dat_i, sp_bio_network, by = "sp")
      } else if(input_format == "matrix"){
        # Merging bioregions of sites with dat
        bioregions_i <- bioregions[, c(1, grep(all_method[i],
                                               colnames(bioregions)))]
        colnames(bioregions_i) <- c("site", "bio_site")
        dat_i <- left_join(sp_df, bioregions_i, by = "site")

        # Assigning bioregions to species
        sp_bio_network <- dat_i %>%
          group_by(sp, bio_site) %>%
          summarise(c = n()) %>%
          filter(row_number(desc(c)) == 1) %>%
          rename(bio_sp = bio_site) %>%
          select(sp, bio_sp)

        dat_i <- left_join(dat_i, sp_bio_network, by = "sp")
      }
    }
    # contribute() function
    scores <- contribute(
      dat = dat_i, sp_col = "sp", site_col = "site",
      bioregion_col = "bio_site", bioregion_sp_col = "bio_sp", ab = NULL)

    # Only columns of interest
    scores <- scores[, c("site", "sp", "bio_site", "bio_sp",
                         "C_site", "C_sp", "z_site", "z_sp", "rho")]
    scores$method <- all_method[i]

    # list element
    contrib_list[[i]] <- scores

    ## interact() ----
    dat_int <- interact(dat = scores, site = "site", sp = "sp",
                        bio_site = "bio_site", plot = TRUE, chord = FALSE)
    # list element
    interact_list[[i]] <- dat_int
  }

  ## Return outputs ----
  if(!is.null(cluster_method)){
    if(!is.null(network_method)){
      return(list(dat, dat_proj, res_cluster, res_network, comp,
                  contrib_list, interact_list))
    } else{
      return(list(dat, dat_proj, res_cluster, comp, contrib_list,
                  interact_list))
    }
  } else{
    return(list(dat, dat_proj, res_network, comp, contrib_list,
                interact_list))
  }
}
