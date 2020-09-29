
contingency <- function(dat, site, sp, ab = NULL, weight = FALSE){

  if(!is.data.frame(dat)){
    stop("dat must be a data.frame with columns sp and site.")
  }

  if(!is.character(site)){
    stop("site must be the column name of dat describing the sites.")
  }

  if(!is.character(sp)){
    stop("sp must be the column name of dat describing the species.")
  }

  if(!is.null(ab) & !is.character(ab)){
    stop("ab must be the column name of dat describing the abundances
         of species.")
  }

  if(is.null(ab) & weight == TRUE){
    warning("Without column abundances, contingency table will only get binary
        values.")
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean.")
  }

  # Rename columns
  colnames(dat)[colnames(dat) == sp] <- "sp"
  colnames(dat)[colnames(dat) == site] <- "site"
  if(!is.null(ab)){
    colnames(dat)[colnames(dat) == ab] <- "ab"
  } else{ # without abundances, species are just present in the sites
    dat$ab <- 1
  }

  # Create contingency table with matrix indexing
  dat$site <- as.factor(dat$site)
  dat$sp <- as.factor(dat$sp)

  mat <- with(dat, {
    out <- matrix(nrow = nlevels(site), ncol = nlevels(sp),
                  dimnames = list(levels(site), levels(sp)))
    out[cbind(site, sp)] <- ab
    out
  })
  # Replace NAs with 0s
  mat[is.na(mat)] <- 0

  # Check for empty rows and columns
  mat <- mat[rowSums(mat) > 0, ]
  mat <- mat[, colSums(mat) > 0]

  # Conversion as binary matrix if weight == FALSE
  if(weight == FALSE){
    mat[mat > 0] <- 1
  }

  return(mat)
}
