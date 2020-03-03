
run_infomap <- function(dat,
                        saving_directory){

  # bioRgeo directory
  Bio_dir <- list.dirs(.libPaths(), recursive = FALSE)
  Bio_dir <- Bio_dir[grep("bioRgeo", Bio_dir)]

  # Change working directory so the file is saved in the proper place
  current_path <- getwd()
  setwd(Bio_dir)

  # Set up the command with required parameters
  if(.Platform$OS.type == "windows"){
    cmd <-
      paste0(Bio_dir,
             "Infomap_win -N 10 --undirected --two-level --map --clu --tree
             --markov-time 0.5 test_uni_inf.txt test")
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
  tp_res <- readLines(paste0(Bio_dir, "/OSLOM/dataset.txt_oslo_files/tp"))

  # Reset previous working directory
  setwd(current_path)

}
