
all_steps <- function(dat, input_format = "tidy", binary = TRUE,
                      site, sp, ab = NULL,
                      metric = "simpson",
                      cluster_method = NULL, network_method = NULL,
                      weight = FALSE,
                      clustering = TRUE, ward_method = "ward.D2",
                      optim_method = "firstSEmax", nstart = 25, B = 50,
                      K.max = 20){
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

    if(is.null(ab) & binary == FALSE){
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

  if(!is.null(cluster_method) &
     !(cluster_method %in% c("kmeans", "meanshift",
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

  if(!is.null(network_method) &
     !(network_method %in% c("greedy", "betweenness", "walktrap", "louvain", "LPAwb",
                             "infomap", "spinglass", "leading_eigen",
                             "label_prop", "netcarto", "oslom", "infomap",
                             "beckett", "quanbimo"))){
    stop("Network algorithm chosen is not available.
     Please chose among the following network algorithms:
         'greedy', 'betweenness', 'walktrap', 'louvain', 'LPAwb',
         'spinglass', 'leading_eigen', 'label_prop', 'netcarto', 'oslom',
         'infomap', 'beckett' or 'quanbimo.'")
  }

  ##
  if(!is.logical(binary)){
    stop("binary must be a boolean.")
  }

  if(!(network_algo %in% c("both", "projected", "bipartite"))){
    stop("network_algo designs the type of algorithm among the following
         choices: projected, bipartite or both.")
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean.")
  }

  if(!(optim_method %in% c("globalmax", "firstmax", "Tibs2001SEmax",
                           "firstSEmax", "globalSEmax."))){
    stop("Chosen gap statistic to determine the optimal number of cluster is
    not available.
     Please chose among the followings:
         globalmax, firstmax, Tibs2001SEmax, firstSEmax or globalSEmax.")
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

  if(form == "contingency" & K.max > nrow(dat)){
    stop("K.max should not be superior to the number of rows of the
         contincenty matrix.")
  }

  ## contingency() ----
  # Create contingency table if needed
  if(input_format == "tidy"){
    sp_mat <- contingency(dat, site = site, sp = sp, ab = ab, binary = binary)
  }

  ## simil() ----
  # Project network with simil()
  dat_proj <- simil(sp_mat, metric = metric, input = "matrix",
                    output = "data frame",
                    site = NULL, sp = NULL, ab = NULL, binary = TRUE)

  dat_proj <- dat_proj[, c("id1", "id2", metric)]

  ## cluster() ----


  ## community() ----

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
