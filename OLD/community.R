
community <- function(dat, algorithm = "greedy", weight = FALSE,
                      input = "matrix", site.col = NULL, species.col = NULL,
                      abundance.col = NULL, N = 10, n_runs = 10, t_param = 0.1,
                      cp_param = 0.5, hr = 0, oslom_id1 = NULL,
                      oslom_id2 = NULL, oslom_proj = NULL){

  ## 1. Controls ----
  if(!(input %in% c("matrix", "data.frame"))){
    stop("dat must be either a contingency matrix with sites as rows and
    species as columns or a long format data.frame with each row being the
    presence of a species in a given site.")
  }

  if(input == "matrix"){
    if(!is.matrix(dat)){
      stop("dat should be a matrix with sites as rows and species as columns.")
    }
    if(is.null(rownames(dat)) | is.null(colnames(dat))){
      stop("dat does not have rownames or colnames.")
    }

  } else if(input == "data.frame"){
    if(!is.data.frame(dat)){
      stop("dat should be a long format data.frame with each row being the
    presence of a species in a given site.")

      if(!is.character(site.col)){
        stop("site.col must be the column name of dat describing the sites")
      }

      if(!is.character(species.col)){
        stop("species.col must be the column name of dat describing the species.")
      }

      if(!is.null(abundance.col) & !is.character(abundance.col)){
        stop("abundance.col must be the column name of dat describing the abundances
         of species")
      }
    }

    if(weight == TRUE & is.null(abundance.col)){
      stop("With weight = TRUE and input = 'data.frame', input data.frame
      should have a column containg the abundances of species at sites.")
    }
  }

  if(!(all(algorithm %in% c("greedy", "betweenness", "walktrap", "louvain",
                            "infomap", "spinglass", "leading_eigen", "label_prop",
                            "netcarto", "oslom", "infomap",
                            "beckett", "quanbimo")))){
    stop("Provided algorithm to detect communities is not available. Please
    choose among the followings: 'greedy', 'betweenness', 'walktrap', 'louvain',
         'spinglass', 'leading_eigen', 'label_prop', 'netcarto', 'oslom',
         'infomap', 'beckett' or 'quanbimo'")
  }

  if(any(algorithm %in% c("beckett", "quanbimo")) & input == "data.frame"){
    warning("beckett and quanbimo algorithms need a site-species matrix as
    an input. We here perform bioRgeo::contingency() function to convert it but
        you can do this step a priori to save time.")
    dat_to_matrix <- contingency(dat, site.col, species.col, abundance.col)
  }

  if(any(algorithm %in% c("infomap")) & input == "matrix"){
    warning("infomap algorithm needs a long-format data.frame as an input.
    We here perform bioRgeo::contingency_to_df() function to convert it but
        you can do this step a priori to save time.")
    dat_to_df <- contingency_to_df(dat, site.col = "site", species.col = "sp",
                                   abundance.col = "ab")
  }

  if(any(algorithm %in% c("netcarto", "quanbimo"))){
    warning("netcarto and quanbimo algorithm is slow. Depending on the size of the
            network, it may take a long time to run.")
  }

  if(!is.logical(weight)){
    stop("'weight' must be a boolean.")
  }

  if(!(abs(N - round(N)) < .Machine$double.eps^0.5)){
    stop("N must be an integer setting the number of outer-most loops to run
         before picking the best solution, for Infomap algorithm only.")
  }

  # 1.1. Controls for OSLOM ----
  if(length(algorithm) > 1 & "oslom" %in% algorithm){
    stop("oslom algorithm needs a different input than the other algorithms.
         You should run community() separately with OSLOM.")
  } else if(length(algorithm) == 1 & algorithm %in% "oslom"){
    # Controls: dat must be a data.frame with three columns containing id1, id2
    # and similarity metric
    if(!(is.data.frame(dat))){
      stop("dat must be a data.frame with three columns containing id1, id2
  and a similarity metric. It is the output of 'bioRgeo::simil()' function")
    }

    if(ncol(dat) < 3){
      stop("dat must be a data.frame with at least three columns containing
      id1, id2 and similarity metric")
    }

    if(!(is.numeric(dat[, 3]))){
      stop("dat must be a data.frame with three columns containing id1, id2
  and similarity metric")
    }

    if(0 %in% dat[, 3]){
      stop("OSLOM needs strictly positive weights to run. Remove the null
         weights from the input data.frame.")
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

    if(!is.character(oslom_id1)){
      stop("'oslom_id1' must be the column name of the first pair of sites.")
    }

    if(!is.character(oslom_id2)){
      stop("'oslom_id2' must be the column name of the second pair of sites.")
    }

    if(!is.character(oslom_proj)){
      stop("'oslom_proj' must be the column name of the similarity metric
           between each pair of sites.")
    }
  }

  # Initiate empty list
  network_list <- list()
  # for-loop over the selected algorithms
  for(i in 1:length(algorithm)){

    ## 2. Infomap ----
    if(algorithm[i] == "infomap"){
      if(any(algorithm %in% c("infomap")) & input == "matrix"){
        sauv <- dat
        dat <- dat_df
      }

      # bioRgeo directory
      Bio_dir <- list.dirs(.libPaths(), recursive = FALSE)
      Bio_dir <- Bio_dir[grep("bioRgeo", Bio_dir)]

      # Add control: only one package should match
      if(length(Bio_dir) > 1){
        stop("Two conflicting versions of bioRgeo seem to coexist.")
      }

      # Convert dat to Pajek format (readable by Infomap)
      if(!is.factor(dat[, species.col])){
        dat[, species.col] <- as.factor(dat[, species.col])
      }
      if(!is.factor(dat[, site.col])){
        dat[, site.col] <- as.factor(dat[, site.col])
      }

      species <- data.frame(species.col = levels(dat[, species.col]),
                            id = 1:length(levels(dat[, species.col])))
      sites <- data.frame(site.col = levels(dat[, site.col]),
                          id = length(levels(dat[, species.col])) +
                            1:length(levels(dat[, site.col])))

      if(!is.null(abundance.col)){
        links <- data.frame(from = species$id[match(dat[, species.col],
                                                    species$species.col)],
                            to = sites$id[match(dat[, site.col],
                                                sites$site.col)],
                            weight = ifelse(rep(length(abundance.col),
                                                nrow(dat)),
                                            dat[, abundance.col],
                                            rep(1, nrow(dat))))
      } else{
        links <- data.frame(
          from = species$id[match(dat[, species.col], species$species.col)],
          to = sites$id[match(dat[, site.col], sites$site.col)],
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
          paste0("./Infomap_lin -N ", N, " --undirected --two-level --map --clu --tree --markov-time 0.5 dat.net output")

        stop("To do")
      } else{
        stop("Windows or Unix distributions only.")
      }
      # Execute the command from R
      system(command = cmd)

      # Reset previous working directory
      setwd(current_path)

      # Import tree file created
      network_lab <- read.table(paste0(Bio_dir, "/INFOMAP/output/dat.tree"),
                                skip = 1, sep = " ", nrows = -1)
      colnames(network_lab) <- c("path", "flow", "node", "node_id")

      # Add column with the first degree module
      network_lab$path <- as.character(network_lab$path)
      network_lab$module <- gsub( ":.*$", "", network_lab$path)

      # Category of node
      if(input == "matrix"){
        network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                                  "site", "sp")
      }else if(input == "data.frame"){
        network_lab$cat <- ifelse(network_lab$node %in% unique(dat$site),
                                  "site", "sp")
      }

      # Remove some columns
      network_lab <- network_lab[, c("node", "module", "cat")]

      # Remove the input dataset
      file.remove(paste0(Bio_dir, "/INFOMAP/dat.net"))

      # Remove output files created (.clu .tree .map)
      file.remove(paste0(Bio_dir, "/INFOMAP/output/",
                         dir(paste0(Bio_dir, "INFOMAP/output/"),
                             pattern = "dat.")))

      # Removing data.frame conversion to save space
      if(any(algorithm %in% c("infomap")) & input == "matrix"){
        dat <- sauv
        rm(sauv)
      }

      ## 3. igraph algorithms ----
    } else if(algorithm[i] %in% c("greedy", "betweenness", "walktrap", "louvain",
                                  "spinglass", "leading_eigen", "label_prop")){
      require(igraph)
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

      } else if(input == "data.frame"){
        if(weight == FALSE){
          network <- graph_from_data_frame(dat[, c(species.col, site.col)],
                                           directed = FALSE)
        } else if(weight == TRUE){
          network <- graph_from_data_frame(dat[, c(species.col, site.col)],
                                           directed = FALSE, weighted = TRUE)
        }
      }

      # Community detection algorithm
      if(algorithm[i] == "greedy"){
        network_mod <- cluster_fast_greedy(network)
      } else if(algorithm[i] == "betweenness"){
        network_mod <- cluster_edge_betweenness(network, modularity = TRUE)
      } else if(algorithm[i] == "walktrap"){
        network_mod <- cluster_walktrap(graph = network)
      } else if(algorithm[i] == "louvain"){
        network_mod <- cluster_louvain(graph = network)
        # } else if(algorithm[i] == "infomap"){
        #   network_mod <- cluster_infomap(graph = network)
      } else if(algorithm[i] == "spinglass"){
        network_mod <- cluster_spinglass(graph = network)
      } else if(algorithm[i] == "leading_eigen"){
        network_mod <- cluster_leading_eigen(graph = network)
      } else if(algorithm[i] == "label_prop"){
        network_mod <- cluster_label_prop(graph = network)
      }

      # Create data.frame with modules
      network_lab <- data.frame(node = names(membership(network_mod)),
                                module = as.character(membership(network_mod)))

      # Add modularity score
      # network_lab$modularity <- modularity(network_mod)

      # Add category of the node
      if(input == "matrix"){
        network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                                  "site", "sp")
      }else if(input == "data.frame"){
        network_lab$cat <- ifelse(network_lab$node %in% unique(dat$site),
                                  "site", "sp")
      }
      ## 4. netcarto ----
    } else if(algorithm[i] == "netcarto"){
      if(input == "data.frame"){
        dat <- contingency(dat, site.col = site.col, species.col = species.col,
                           abundance.col = abundance.col)
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
      # network_lab$modularity <- tmp[[2]]

      # Add category of the node
      if(input == "matrix"){
        network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                                  "site.col", "species.col")
      }else if(input == "data.frame"){
        network_lab$cat <- ifelse(network_lab$node %in% dat[, site.col],
                                  "site.col", "species.col")
      }
      ## 5. OSLOM ----
    } else if(algorithm[i] == "oslom"){
      warning("With OSLOM algorithm, the input must be a projected network.")

      # Renaming columns (id1, id2 must be present, no more than 3 columns)
      oslom_dat <- dat[, c(oslom_id1, oslom_id2, oslom_proj)]
      colnames(oslom_dat) <- c("id1", "id2", "oslom_proj")

      # bioRgeo directory
      Bio_dir <- list.dirs(.libPaths(), recursive = FALSE)
      Bio_dir <- Bio_dir[grep("bioRgeo", Bio_dir)]

      # Add control: only one package should match
      if(length(Bio_dir) > 1){
        stop("Two conflicting versions of bioRgeo seem to coexist.")
      }

      # Reformating input data.frame
      oslom_dat <- dat[, c(oslom_id1, oslom_id2, oslom_proj)]
      # Renaming columns (id1, id2 must be present, no more than 3 columns)
      colnames(oslom_dat) <- c("id1", "id2", "oslom_proj")

      # All the columns of oslom_dat have to be numeric
      if(!is.numeric(oslom_dat$id1)){
        oslom_dat$id1 <- as.numeric(as.factor(oslom_dat$id1))
      }
      if(!is.numeric(oslom_dat$id2)){
        oslom_dat$id2 <- as.numeric(as.factor(oslom_dat$id2))
      }

      # Save input dataset as a .txt file into OSLOM folder
      write.table(oslom_dat,
                  paste0(Bio_dir, "/OSLOM/dataset.txt"), row.names = FALSE)

      # Change working directory so the file is saved in the proper place
      current_path <- getwd()
      setwd(Bio_dir)

      # Set up the command with required parameters
      if(.Platform$OS.type == "windows"){
        cmd <-
          paste0(Bio_dir, "/OSLOM/oslom_undir_win.exe -f OSLOM/dataset.txt -w",
                 " -r ", n_runs, " -t ", t_param, " -cp ", cp_param, " -hr ", hr)
      } else if(.Platform$OS.type == "unix"){
        cmd <-
          paste0(Bio_dir, "/OSLOM/oslom_undir_lin -f OSLOM/dataset.txt -w",
                 " -r ", n_runs, " -t ", t_param, " -cp ", cp_param, " -hr ", hr)
      } else{
        stop("Windows or Unix distributions only.")
      }

      # Execute the command from R
      system(command = cmd)

      # Control: if the command line did not work, previous working directory reset
      if(!("tp" %in% list.files(paste0(Bio_dir,
                                       "/OSLOM/dataset.txt_oslo_files")))){
        setwd(current_path)
        stop("Command line was wrongly implemented. OSLOM did not run.")
      }

      # Import tp file created
      tp_res <- readLines(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/tp"))

      # Reset previous working directory
      setwd(current_path)

      # Remove .oslo_files created and the dataset
      file.remove(paste0(Bio_dir, "/OSLOM/dataset.txt"))
      file.remove(paste0(Bio_dir, "/time_seed.dat"))

      # Remove all filed in .oslo_files folder
      file.remove(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/",
                         dir(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/"),
                             pattern = "net")))
      file.remove(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/",
                         dir(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/"),
                             pattern = "partitions_level")))
      file.remove(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/",
                         dir(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/"),
                             pattern = "statistics_level")))
      file.remove(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/",
                         dir(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/"),
                             pattern = "short_tp")))
      file.remove(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/",
                         dir(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/"),
                             pattern = "tp")))

      # OSLOM output
      # Convert tp_res into list
      bioregion_list <- list()
      length(bioregion_list) <- length(tp_res)/2
      for(k in 1:length(tp_res)){
        if((k/2-trunc(k/2)) == 0){ # loop over module pixels
          bioregion_list[[(k/2)]] <-
            as.numeric(as.matrix(strsplit(tp_res[k], split = " ")[[1]]))
        }
      }

      # Convert bioregion_list into data.frame
      names(bioregion_list) <- as.character(seq_along(bioregion_list))
      network_lab <- data.frame(
        id_oslom = unlist(bioregion_list),
        module = rep(names(bioregion_list), sapply(bioregion_list, length)))

      # Column with the names of the sites
      network_lab$site <- unique(
        c(levels(as.factor(sp_proj$id1)),
          levels(as.factor(sp_proj$id2))))[network_lab$id_oslom]

      # If some doublons are present => warning
      if(length(bioregion_list) > 1 &
         length(unique(table(network_lab$id))) > 1){
        doublon <- table(network_lab$id)[table(network_lab$id) > 1]
        warning(paste0(length(doublon), " sites have several bioregions assigned."))
      }

      network_lab <- network_lab[, c("site", "module")]

      ## 6. Beckett/Dormann & Strauss ----
    } else if(algorithm[i] %in% c("beckett", "quanbimo")){
      require(bipartite)
      if(algorithm[i] == "beckett"){
        algorithm2 <- "Beckett"
      } else if(algorithm[i] == "quanbimo"){
        algorithm2 <- "DormannnStrauss"
      }

      if(any(algorithm %in% c("beckett", "quanbimo")) & input == "data.frame"){
        network_lab <- computeModules(
          dat_to_matrix, method = algorithm2, deep = FALSE,
          deleteOriginalFiles = TRUE, steps = 1000000, tolerance = 1e-10,
          experimental = FALSE, forceLPA = FALSE)

        site_names <- rownames(dat_to_matrix)
        sp_names <- colnames(dat_to_matrix)
      } else{
        network_lab <- computeModules(
          dat, method = algorithm2, deep = FALSE,
          deleteOriginalFiles = TRUE, steps = 1000000, tolerance = 1e-10,
          experimental = FALSE, forceLPA = FALSE)

        site_names <- rownames(dat)
        sp_names <- colnames(dat)
      }

      # First row and two first columns ignored (see ?`moduleWeb-class`)
      network_lab <- network_lab@modules[-1, -c(1, 2)]
      # Modules in rows; columns have first the sites and then the species
      colnames(network_lab) <- c(site_names, sp_names)
      # numbers just indicate in which module a node is
      network_lab[network_lab > 0] <- 1
      # conversion to data.frame
      network_lab <- reshape2::melt(network_lab)
      network_lab <- network_lab[which(network_lab$value > 0), ]
      network_lab <- data.frame(node = network_lab$Var2,
                                module = network_lab$Var1)
      network_lab$cat <- ifelse(network_lab$node %in% site_names, "site", "sp")
    }

    # Convert module column into character
    network_lab$module <- as.character(network_lab$module)

    # Changing column names of network_lab: paste(algorithm[i], method)
    colnames(network_lab)[colnames(network_lab) == "module"] <-
      paste0("module_", algorithm[i])

    # Complete list
    network_list[[i]] <- network_lab
  }

  # Combine all list elements into one data.frame
  network_list <- Reduce(merge, network_list)
  return(network_list)
}
