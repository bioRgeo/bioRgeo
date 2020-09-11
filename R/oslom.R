
run_oslom <- function(dat, n_runs = 10, t_param = 0.1, cp_param = 0.5, hr = 0){
  ## 1. Controls ----
  # Controls: dat must be a data.frame with three columns containing id1, id2
  # and similarity metric
  if(!(is.data.frame(dat))){
    stop("dat must be a data.frame with three columns containing id1, id2
  and a similarity metric. It is the output of 'bioRgeo::simil()' function")
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

  ## 2. Running OSLOM ----
  # bioRgeo directory
  Bio_dir <- list.dirs(.libPaths(), recursive = FALSE)
  Bio_dir <- Bio_dir[grep("bioRgeo", Bio_dir)]

  # Add control: only one package should match
  if(length(Bio_dir) > 1){
    stop("Two conflicting versions of bioRgeo seem to coexist.")
  }

  # Save input dataset as a .txt file into OSLOM folder
  write.table(dat, paste0(Bio_dir, "/OSLOM/dataset.txt"), row.names = FALSE)

  # Change working directory so the file is saved in the proper place
  current_path <- getwd()
  setwd(Bio_dir)

  # Set up the command with required parameters
  if(.Platform$OS.type == "windows"){
    cmd <-
      paste0(Bio_dir, "/OSLOM/oslom_undir_win.exe -f OSLOM/dataset.txt -w",
             " -r ", n_runs, " -t ", t_param, " -cp ", cp_param, " -hr ", hr)
  } else if(.Platform$OS.type == "unix"){
    stop("To do")
  } else{
    stop("Windows or Unix distributions only.")
  }

  # Execute the command from R
  system(command = cmd)

  # Control: if the command line did not work, previous working directory reset
  if(!("tp" %in% list.files(paste0(Bio_dir,
                                   "/OSLOM/dataset.txt_oslo_files")))){
    stop("Command line was wrongly implemented. OSLOM did not run.")
    setwd(current_path)
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

  ## 3. OSLOM output ----
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
  bioregion_oslom <- data.frame(
    id_oslom = unlist(bioregion_list),
    bioregion = rep(names(bioregion_list), sapply(bioregion_list, length)))

  # Column with the names of the sites
  bioregion_oslom$site <- unique(
    c(levels(as.factor(sp_proj$id1)),
      levels(as.factor(sp_proj$id2))))[bioregion_oslom$id_oslom]

  # If some doublons are present => warning
  if(length(bioregion_list) > 1 &
     length(unique(table(bioregion_oslom$id))) > 1){
    doublon <- table(bioregion_oslom$id)[table(bioregion_oslom$id) > 1]
    warning(paste0(length(doublon), " sites have several bioregions assigned."))
  }

  return(bioregion_oslom[, c("site", "bioregion")])
}
