
contingency_to_df <- function(dat, site.col = "site", species.col = "species",
                              abundance.col = "abundances",
                              remove_absent_species = TRUE){
  if(!is.matrix(dat)){
    stop("Contingency table should be a matrix with sites as rows and
         species as columns.")
  }

  if(!is.character(site.col)){
    stop("site.col must be the column name describing the rows of the contingency
         matrix.")
  }

  if(!is.character(species.col)){
    stop("species.col must be the column name describing the columns of the
    contingency matrix.")
  }

  if(!is.character(abundance.col)){
    stop("abundance.col must be the column name describing the values of the
    contingency matrix.")
  }

  if(!is.logical(remove_absent_species)){
    stop("remove_absent_species must be a boolean determining whether absent
    species from the contingency matrix have to be removed from the output.")
  }

  # Conversion as data.frame
  dat_df <- reshape2::melt(dat)
  colnames(dat_df) <- c(site.col, species.col, abundance.col)

  if(remove_absent_species == TRUE){
    dat_df <- dat_df[which(dat_df$ab > 0), ]
  }

  return(dat_df)
}
