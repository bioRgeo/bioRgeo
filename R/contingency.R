
contingency <- function(dat, site.col, species.col, abundance.col = NULL){

  # Control 1
  if(!is.data.frame(dat)){
    stop("dat must be a data.frame with columns species and sites.")
  }

  if(!is.character(site.col)){
    stop("site.col must be the column name of dat describing the sites.")
  }

  if(!is.character(species.col)){
    stop("species.col must be the column name of dat describing the species.")
  }

  if(!is.null(abundance.col) & !is.character(abundance.col)){
    stop("abundance.col must be the column name of dat describing the
    abundances of species.")
  }

  # Control 2
  if(sum(colnames(dat) == site.col)==0){
    stop(paste0("could not find the column ", species.col))
  }

  if(sum(colnames(dat) == species.col)==0){
    stop(paste0("could not find the column ", site.col))
  }

  if(!is.null(abundance.col) & sum(colnames(dat) == abundance.col)==0){
    stop(paste0("could not find the column ", abundance.col))
  }

  # Rename columns
  colnames(dat)[colnames(dat) == species.col] <- "species.col"
  colnames(dat)[colnames(dat) == site.col] <- "site.col"
  if(!is.null(abundance.col)){
    colnames(dat)[colnames(dat) == abundance.col] <- "abundance.col"
  } else{ # without abundances, species are just present in the sites
    dat$abundance.col <- 1
  }

  # Create contingency table with matrix indexing
  dat$site.col <- as.factor(dat$site.col)
  dat$species.col <- as.factor(dat$species.col)

  mat <- with(dat, {
    out <- matrix(nrow = nlevels(site.col), ncol = nlevels(species.col),
                  dimnames = list(levels(site.col), levels(species.col)))
    out[cbind(site.col, species.col)] <- abundance.col
    out
  })
  # Replace NAs with 0s
  mat[is.na(mat)] <- 0

  # Check for empty rows and columns
  mat <- mat[rowSums(mat) > 0, ]
  mat <- mat[, colSums(mat) > 0]

  # Conversion as binary matrix if is.null(abundance.col)
  if(is.null(abundance.col)){
    mat[mat > 0] <- 1
  }

  return(mat)
}
