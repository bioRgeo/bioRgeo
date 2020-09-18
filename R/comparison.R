
comparison <- function(dat, site, bio_col, thres = 10, output = "both"){

  ## 1. Controls ----
  if(!is.data.frame(dat)){
    stop("Input must be a data.frame with a column containing the sites and
         columns containing different partitions.")
  }

  if(output %in% c("common", "both")){
    if(!(site %in% colnames(dat))){
      stop("dat must contains a column 'site' with all the sites.")
    }
  }

  if(!is.numeric(bio_col) | max(bio_col) > ncol(dat)){
    stop("bio_col must be a numeric indicating which columns of dat contains
         the different bioregions.")
  }

  if(!is.numeric(thres)){
    stop("thres must be a numeric determining the minimum number of pixels
         per pair.")
  }

  if(!(output %in% c("percentage", "common", "both"))){
    stop("output must be a character string indicating whether the output is
         is a square matrix with the number of times each pair of sites have
         been classed together ('percentage'), a data frame with the sites all
         the time classed together ('common') or both objects ('both').")
  }

  ## 2. Groups of sites ----
  # Get matrix with number of bioregions in common per pair of pixels
  list_mat_pair <- list()
  for(i in 1:length(bio_col)){
    # matrix with sites in rows & columns
    # TRUE when both sites have same bioregion
    list_mat_pair[[i]] <-
      matrix(outer(dat[, bio_col[i]], dat[, bio_col[i]], "=="),
             nrow = nrow(dat), dimnames = list(dat[, site], dat[, site]))
  }
  # Sum all matrices: how many times each pair of sites has same bioregion
  list_mat_pair <- Reduce('+', list_mat_pair)

  # Remove lower part and diagonal of the matrix
  diag(list_mat_pair) <- NA
  list_mat_pair[lower.tri(list_mat_pair)] <- NA

  if(output == "percentage"){
    return(list_mat_pair)
  } else if(output %in% c("common", "both")){

    # percentage instead of sum
    list_mat_pair <- 100*list_mat_pair/length(bio_col)

    # Conversion to data frame
    list_df_pair <- reshape2::melt(list_mat_pair)
    colnames(list_df_pair) <- c("id1", "id2", "perc")

    # Remove NAs
    list_df_pair <- list_df_pair[complete.cases(list_df_pair), ]

    # Extract pixels together 100% of the time
    all100 <- list_df_pair[which(list_df_pair$perc == 100), ]

    # Add bioregion (whatever algorithm since all pixels are grouped together)
    all100 <- left_join(all100,
                        dat[, c(site, colnames(dat)[bio_col[1]])],
                        by = c("id1" = site))
    colnames(all100) <- c("id1", "id2", "perc", "uni")

    # Rename levels of common bioregion
    all100$uni <- as.factor(all100$uni)
    levels(all100$uni) <- as.character(seq(1:length(levels(all100$uni))))

    # Remove pairs of pixels under a given threshold of pixels
    thresh <- names(table(all100$uni)[table(all100$uni) > thres])

    all100 <- all100[which(all100$uni %in% thresh), ]

    if(output == "common"){
      return(all100)
    } else if(output == "both"){
      return(list(list_mat_pair, all100))
    }
  }
}
