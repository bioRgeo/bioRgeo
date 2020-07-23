
algo_bipartite <- function(dat, algo = "greedy", weight = FALSE,
                           input = "matrix", site = NULL, sp = NULL,
                           ab = NULL, saving_directory, N = 10){

  ## 1. Controls ----
  if(!(input %in% c("matrix", "data frame"))){
    stop("dat must be either a contingency matrix with sites as rows and
    species as columns or a long format data frame with each row being the
    presence of a species in a given site.")
  }

  if(input == "matrix"){
    if(!is.matrix(dat)){
      stop("dat should be a matrix with sites as rows and species as columns.")
    }
  } else if(input == "data frame"){
    if(!is.data.frame(dat)){
      stop("dat should be a long format data frame with each row being the
    presence of a species in a given site.")

      if(!is.character(site)){
        stop("site must be the column name of dat describing the sites")
      }

      if(!is.character(sp)){
        stop("sp must be the column name of dat describing the species.")
      }

      if(!is.null(ab) & !is.character(ab)){
        stop("ab must be the column name of dat describing the abundances
         of species")
      }
    }

    if(weight == TRUE & is.null(ab)){
      stop("With weight = TRUE and input = 'data frame', input data frame
      should have a column containg the abundances of species at sites.")
    }
  }

  if(!(algo %in% c("greedy", "girvan", "walktrap", "louvain", "LPAwb",
                   "infomap", "spinglass", "leading_eigen", "label_prop",
                   "optimal", "netcarto"))){
    stop("Provided algorithm to compute modularity is not available. Please
    choose among the followings: greedy, girvan, walktrap, louvain, LPAwb,
         infomap, spinglass, leading_eigen, label_prop, optimal or
         netcarto.")
  }

  if(algo %in% c("optimal", "netcarto")){
    warning("The chosen algorithm may take time to run.")
  }

  if(!is.logical(weight)){
    stop("'weight' must be a boolean.")
  }

  if(!(is.character(saving_directory)) |
     !dir.exists(saving_directory)){
    stop("saving_directory must be a path where the Infomap file containing
         the bioregions identified will be saved.")
  }

  if(!(abs(N - round(N)) < .Machine$double.eps^0.5)){
    stop("N must be an integer setting the number of outer-most loops to run
         before picking the best solution.")
  }

  require(igraph)
  # source("scripts/Beckett_LPAwb_function.R") => store in R function

  ## 2. Infomap ----
  if(algo == "infomap"){
    # bioRgeo directory
    Bio_dir <- list.dirs(.libPaths(), recursive = FALSE)
    Bio_dir <- Bio_dir[grep("bioRgeo", Bio_dir)]

    # Add control: only one package should match
    if(length(Bio_dir) > 1){
      stop("Two conflicting versions of bioRgeo seem to coexist.")
    }

    # Convert dat to Pajek format (readable by Infomap)
    if(!is.factor(dat[, sp])){
      dat[, sp] <- as.factor(dat[, sp])
    }
    if(!is.factor(dat[, site])){
      dat[, site] <- as.factor(dat[, site])
    }

    species <- data.frame(sp = levels(dat[, sp]),
                          id = 1:length(levels(dat[, sp])))
    sites <- data.frame(site = levels(dat[, site]),
                        id = length(levels(dat[, sp])) +
                          1:length(levels(dat[, site])))

    if(!is.null(ab)){
      links <- data.frame(from = species$id[match(dat[, sp], species$sp)],
                          to = sites$id[match(dat[, site], sites$site)],
                          weight = ifelse(rep(length(ab), nrow(dat)),
                                          dat[, ab],
                                          rep(1, nrow(dat))))
    } else{
      links <- data.frame(from = species$id[match(dat[, sp], species$sp)],
                          to = sites$id[match(dat[, site], sites$site)],
                          weight = rep(1, nrow(dat)))
    }

    # Write .net file into INFOMAP folder
    cat(paste0("*Vertices ", max(sites$id), "\n",
               paste0(species$id, ' "', species$sp, '"', collapse = "\n"), "\n",
               paste0(sites$id, ' "', sites$site, '"', collapse = "\n"),
               "\n", "*Edges\n",
               paste(apply(links, 1, paste, collapse = " "), collapse = "\n")),
        file = paste0(Bio_dir, "/INFOMAP/dat.net"))

    # Change working directory so the file is saved in the proper place
    current_path <- getwd()
    setwd(paste0(Bio_dir, "/INFOMAP"))

    # Set up the command with required parameters
    if(.Platform$OS.type == "windows"){
      cmd <- paste0("Infomap_win -N ", N, " --undirected --two-level --map --clu --tree --markov-time 1 dat.net output")

      # "  Infomap_win -N 10 --undirected --two-level --skip-adjust-bipartite-flow
      # --map --clu -i'bipartite' test_bip_inf.txt test"

    } else if(.Platform$OS.type == "unix"){
      cmd <-
        paste0(Bio_dir,
               "./Infomap_lin -N 10 --undirected --two-level --map --clu --tree --markov-time 0.5 dataset.txt test")

      stop("To do")
    } else{
      stop("Windows or Unix distributions only.")
    }
    # Execute the command from R
    system(command = cmd)

    # Reset previous working directory
    setwd(current_path)

    # Import tree file created
    tree_infomap <- read.table(paste0(Bio_dir, "/INFOMAP/output/dat.tree"),
                               skip = 1, sep = " ", nrows = -1)
    colnames(tree_infomap) <- c("path", "flow", "name", "node_id")
    tree_infomap$path <- as.character(tree_infomap$path)

    # Add column with the first degree module
    tree_infomap$module <- gsub( ":.*$", "", tree_infomap$path)

    # Category of node
    tree_infomap$cat <- ifelse(tree_infomap$name %in% unique(dat$site),
                               "site", "sp")

    # Merging modules of the site with the input
    infomap_site <- tree_infomap[which(tree_infomap$cat == site),
                                 c("name", "module")]
    colnames(infomap_site) <- c(site, "module")
    dat <- dplyr::left_join(dat, infomap_site, by = site)

    # Merging modules of the species with the input
    infomap_sp <- tree_infomap[which(tree_infomap$cat == sp),
                               c("name", "module")]
    colnames(infomap_sp) <- c(sp, "module_sp")
    dat <- dplyr::left_join(dat, infomap_sp, by = sp)

    # Saving .tp file with bioregions into chosen saving_directory
    saveRDS(dat, file = paste0(saving_directory, "/dat_infomap.rds"))

    # Remove the input dataset
    file.remove(paste0(Bio_dir, "/INFOMAP/dat.net"))

    # Remove output files created (.clu .tree .map)
    file.remove(paste0(Bio_dir, "/INFOMAP/output/",
                       dir(paste0(Bio_dir, "INFOMAP/output/"),
                           pattern = "dat.")))

    network_lab <- infomap_site

    ## 3. LPAwb ----
  } else if(algo == "LPAwb"){ # Beckett modularity
    dat <- as.matrix(dat)
    # Find labels and weighted modularity using LPAwb+
    network_mod <- bioRgeo::LPA_wb_plus(dat)

    # Conversion into data.frame with node, category and module
    network_lab <- data.frame(node = c(rownames(dat), colnames(dat)),
                              module = c(as.character(network_mod[[1]]),
                                         as.character(network_mod[[2]])),
                              modularity = as.character(network_mod[[3]]))
    # Add category of the node
    network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                              "site", "sp")

    ## 4. igraph algorithms ----
  } else if(algo %in% c("greedy", "girvan", "walktrap", "louvain",
                        "spinglass", "leading_eigen", "label_prop",
                        "optimal")){
    # https://stats.stackexchange.com/questions/209086/community-detection-and-modularity
    # https://www.sixhat.net/finding-communities-in-networks-with-r-and-igraph.html

    if(input == "matrix"){
      # Convert matrix into square matrix (first pixels and then species)
      rownames_mat <- rownames(dat)
      dat_sq <- rbind(
        cbind(array(0, c(nrow(dat), nrow(dat))), as.matrix(dat)),
        cbind(as.matrix(t(dat)), array(0, c(ncol(dat), ncol(dat)))))

      # Add pixel names to square matrix
      rownames(dat_sq)[1:length(rownames_mat)] <- rownames_mat
      colnames(dat_sq)[1:length(rownames_mat)] <- rownames_mat

      if(weight == FALSE){
        # Tranforming site_sp matrix into binary matrix
        dat_sq[dat_sq > 0] <- 1
        network <- graph_from_adjacency_matrix(
          dat_sq, mode = "undirected", add.rownames = TRUE, weighted = NULL)
      } else if(weight == TRUE){
        network <- graph_from_adjacency_matrix(
          dat_sq, mode = "undirected", add.rownames = TRUE, weighted = TRUE)
      }

    }else if(input == "data frame"){
      if(weight == FALSE){
        network <- graph_from_data_frame(dat[, c(sp, site)], directed = FALSE)
      } else if(weight == TRUE){
        network <- graph_from_data_frame(dat[, c(sp, site)], directed = FALSE,
                                         weighted = TRUE)
      }
    }

    # Modularity algorithm
    if(algo == "greedy"){
      network_mod <- cluster_fast_greedy(network)
    } else if(algo == "girvan"){
      network_mod <- cluster_edge_betweenness(network, modularity = TRUE)
    } else if(algo == "walktrap"){
      network_mod <- cluster_walktrap(graph = network)
    } else if(algo == "louvain"){
      network_mod <- cluster_louvain(graph = network)
    } else if(algo == "infomap"){
      network_mod <- cluster_infomap(graph = network)
    } else if(algo == "spinglass"){
      network_mod <- cluster_spinglass(graph = network)
    } else if(algo == "leading_eigen"){
      network_mod <- cluster_leading_eigen(graph = network)
    } else if(algo == "label_prop"){
      network_mod <- cluster_label_prop(graph = network)
    } else if(algo == "optimal"){
      network_mod <- cluster_optimal(graph = network)
    }

    # Create data.frame with modules
    network_lab <- data.frame(node = names(membership(network_mod)),
                              module = as.character(membership(network_mod)))

    # Add modularity score
    network_lab$modularity <- modularity(network_mod)

    # Add category of the node
    if(input == "matrix"){
      network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                                "site", "sp")
    }else if(input == "data frame"){
      network_lab$cat <- ifelse(network_lab$node %in% dat[, site],
                                "site", "sp")
    }
    ## 5. netcarto ----
  } else if(algo == "netcarto"){
    if(input == "data frame"){
      dat <- contingency(dat, site = site, sp = sp, ab = ab)
    }
    # Convert matrix into square matrix (first pixels and then species)
    rownames_mat <- rownames(dat)
    dat_sq <- rbind(
      cbind(array(0, c(nrow(dat), nrow(dat))), as.matrix(dat)),
      cbind(as.matrix(t(dat)), array(0, c(ncol(dat), ncol(dat)))))

    # Add pixel names to square matrix
    rownames(dat_sq)[1:length(rownames_mat)] <- rownames_mat
    colnames(dat_sq)[1:length(rownames_mat)] <- rownames_mat
    tmp <- rnetcarto::netcarto(web = dat_sq)

    # Create data.frame with modules
    network_lab <- data.frame(node = tmp[[1]][, "name"],
                              module = tmp[[1]][, "module"])

    # Add modularity score
    network_lab$modularity <- tmp[[2]]

    # Add category of the node
    if(input == "matrix"){
      network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                                "site", "sp")
    }else if(input == "data frame"){
      network_lab$cat <- ifelse(network_lab$node %in% dat[, site],
                                "site", "sp")
    }

  }

  return(network_lab)

  # Also try with bipartite package
  # computeModules(web, method = "Beckett", deep = FALSE,
  #                deleteOriginalFiles = TRUE, steps = 1000000, tolerance = 1e-10,
  #                experimental = FALSE, forceLPA = FALSE)

}
