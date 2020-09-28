
interact <- function(dat, site, sp, bio_site, plot = FALSE, chord = FALSE){

  # Controls ------------------------------------------------------------------
  if(!is.data.frame(dat)){
    stop("dat must be a data.frame with four columns containing the sites,
         the species occurrences in each site and the bioregion of each site.")
  }

  if(!is.character(site)){
    stop("site must be the column name of dat describing the sites.")
  }

  if(!is.character(sp)){
    stop("sp must be the column name of dat describing the species occurrences
         within sites.")
  }

  if(!is.character(bio_site)){
    stop("bio_site must be the column name of dat describing the bioregions
         of each site.")
  }

  if(!is.logical(plot)){
    stop("plot argument should be a boolean determining whether the heatmap of
         interactions between bioregions should be saved.")
  }

  if(!is.logical(chord)){
    stop("chord argument should be a boolean determining whether the chord diagram of
         interactions between bioregions should be saved.")
  }

  # Sorensen ------------------------------------------------------------------
  beta_bio <- c()

  # data.frame with each pair of bioregions
  bio_comb <- data.frame(expand.grid(sort(unique(dat[, bio_site])),
                                     sort(unique(dat[, bio_site]))))
  bio_comb <- data.frame(bioregion1 = as.character(bio_comb$Var1),
                         bioregion2 = as.character(bio_comb$Var2))

  # loop over the pairs of bioregions
  for(i in 1:nrow(bio_comb)){
    bio_i <- unique(dat[which(dat[, bio_site] == bio_comb[i, "bioregion1"]),
                        sp])
    bio_j <- unique(dat[which(dat[, bio_site] == bio_comb[i, "bioregion2"]),
                        sp])

    aij <- table(bio_i %in% bio_j)
    aij <- as.numeric(aij[names(aij) == "TRUE"]) # common sites
    if(length(aij) == 0){
      common_ij <- 0
    } else{
      common_ij <- aij
    }
    # Binding results
    beta_bio <- rbind(
      beta_bio,
      data.frame(
        bio_i = as.character(bio_comb[i, 1]),
        bio_j = as.character(bio_comb[i, 2]),
        beta = 2*common_ij / sum(length(bio_i) + length(bio_j))))
  }

  # Conversion to square matrix to delete duplicated combinations
  beta_bio_mat <- bioRgeo::contingency(beta_bio, site = "bio_i", sp = "bio_j",
                                       ab = "beta", binary = FALSE)
  # removal of upper diagonal
  beta_bio_mat[upper.tri(beta_bio_mat)] <- NA
  # conversion to data.frame
  beta_bio2 <- reshape2::melt(beta_bio_mat)
  colnames(beta_bio2) <- c("Bioregion_1", "Bioregion_2", "beta")

  if(plot == TRUE){ # Heatmap
    require(ggplot2)
    res_plot <-  ggplot(beta_bio2, aes(Bioregion_1, Bioregion_2)) +
      geom_tile(aes(fill = beta), color = "white") +
      geom_text(aes(label = round(beta, 2)), color = "black", size = 6) +
      scale_fill_gradient("Sorensen",
                          #low = "#023858", high = "#a6bddb",
                          low = "grey90", high = "grey20",
                          na.value = "white") +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA),
            legend.justification = c(1, 0),
            legend.position = c(0.4, 0.85),
            legend.direction = "horizontal")
  }

  if(chord == TRUE){ # Chord diagram
    require(circlize)
    circ_beta <- beta_bio2[complete.cases(beta_bio2), ]

    # color gradient
    if(length(unique(dat[, bio_site])) <= 2){
      grid_col <- c("firebrick3", "dodgerblue")
      names(grid_col) <- sort(unique(dat[, bio_site]))
    }else if(length(unique(dat[, bio_site])) <= 9){
      require(RColorBrewer)
      grid_col <- brewer.pal(length(unique(dat[, bio_site])),
                             name = "Paired")
      names(grid_col) <- sort(unique(dat[, bio_site]))
    } else{
      grid_ramp <- colorRampPalette(c("grey", "black"))
      grid_col <- grid_ramp(length(unique(dat[, bio_site])))
      names(grid_col) <- sort(unique(dat[, bio_site]))
    }

    chordDiagram(circ_beta, #[which(circ_beta$Bioregion_1 != circ_beta$Bioregion_2), ],
                 transparency = 0.1,
                 col = grid_col[as.character(circ_beta$Bioregion_1)],
                 grid.col = grid_col,
                 annotationTrack = c("name", "grid"), directional = TRUE,
                 link.border = "black", grid.border = "black")
    title("Sorensen similarity between bioregions")
  }

  # Output
  if(plot == FALSE){
    return(beta_bio = beta_bio2)
  } else if(plot == TRUE){
    return(list(beta_bio = beta_bio2, res_plot))
  }
}
