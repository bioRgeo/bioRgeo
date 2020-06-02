
run_infomap <- function(dat, site = NULL, sp = NULL, ab = NULL,
                        saving_directory, N = 10){

  if(!is.data.frame(dat)){
    stop("dat should be a long format data frame with each row being the
    presence of a species in a given site.")
  }

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

  if(!(is.character(saving_directory)) |
     !dir.exists(saving_directory)){
    stop("saving_directory must be a path where the OSLOM .tp file containing
         the bioregions identified will be saved.")
  }

  if(!(abs(N - round(N)) < .Machine$double.eps^0.5)){
    stop("N must be an integer setting the number of outer-most loops to run
         before picking the best solution.")
  }

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
}
