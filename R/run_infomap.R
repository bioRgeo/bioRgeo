
run_infomap <- function(dat, weight = FALSE, site = NULL, sp = NULL,
                        ab = NULL, saving_directory){

  if(!is.data.frame(dat)){
    stop("dat should be a long format data frame with each row being the
    presence of a species in a given site.")
  }

  if(weight == TRUE & is.null(ab)){
    stop("With weight = TRUE and input = 'data frame', input data frame
      should have a column containg the abundances of species at sites.")
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

  # bioRgeo directory
  Bio_dir <- list.dirs(.libPaths(), recursive = FALSE)
  Bio_dir <- Bio_dir[grep("bioRgeo", Bio_dir)]

  # Convert dat to Pajek format (readable by Infomap)
  # if(!is.factor(dat[, sp])){
  #   dat[, sp] <- as.factor(dat[, sp])
  # }
  # if(!is.factor(dat[, site])){
  #   dat[, site] <- as.factor(dat[, site])
  # }
  #
  # species <- data.frame(sp = levels(dat[, sp]),
  #                       id = 1:length(levels(dat[, sp])))
  # sites <- data.frame(site = levels(dat[, site]),
  #                     id = length(levels(dat[, sp])) +
  #                       1:length(levels(dat[, site])))
  #
  # links <- data.frame(from = species$id[match(dat[, sp], species$sp)],
  #                     to = sites$id[match(dat[, site], sites$site)],
  #                     weight = ifelse(rep(length(ab), nrow(dat)),
  #                                     dat[, ab],
  #                                     rep(1, nrow(dat))))
  #
  # cat(paste0("*Vertices ", max(sites$id), "\n",
  #            paste0(species$id, ' "', species$sp, '"', collapse = "\n"), "\n",
  #            paste0(sites$id, ' "', sites$site, '"', collapse = "\n"),
  #            "\n", "*Edges\n",
  #            paste(apply(links, 1, paste, collapse = " "), collapse = "\n")),
  #     file = "dat_Pajek.net")

  # Save input dataset as a .txt file into INFOMAP folder
  write.table(dat, paste0(Bio_dir, "/INFOMAP/dataset.txt"), row.names = FALSE)

  # Change working directory so the file is saved in the proper place
  current_path <- getwd()
  setwd(Bio_dir)

  # Set up the command with required parameters
  if(.Platform$OS.type == "windows"){
    cmd <-
      paste0(Bio_dir,
             "Infomap_win -N 10 --undirected --two-level --map --clu --tree
             --markov-time 0.5 dat_Pajek.net test")
    # "  Infomap_win -N 10 --undirected --two-level --skip-adjust-bipartite-flow
    # --map --clu -i'bipartite' test_bip_inf.txt test"

  } else if(.Platform$OS.type == "unix"){
    stop("To do")
  } else{
    stop("Windows or Unix distributions only.")
  }
  # Execute the command from R
  system(command = cmd)

  # Import tp file created


  # Reset previous working directory
  setwd(current_path)

}
