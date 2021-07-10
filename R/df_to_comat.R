#' Create a co-occurence matrix from a data.frame
#'
#' This function creates a co-occurence matrix with sites as rows and species as columns
#' from a two- or three-columns data.frame where
#' each row represents the presence of a species (second column)
#' in a site (first column) and an optional third column indicating the abundance
#' of the species in this site (if weight=TRUE).
#'
#' @param df a two- or three-columns data.frame where each row represents the presence of
#' a species (second column) in a site (first column) and an optional third column indicating
#' the abundance of the species in this site
#' @param weight a boolean indicating if the abundance should be considered
#' @export
#' @return A data.frame  where each row represents the presence of a species (second column)
#' in a site (first column) and an optional third column indicating the abundance
#' of the species in this site
#' @author
#' Pierre Denelle \email{pierre.denelle@gmail.com}
#' Maxime Lenormand \email{maxime.lenormand@inrae.fr}
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' df <- data.frame(Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#' Species = c("a", "b", "a", "c", "d", "b", "d"),
#' Weight = c(10, 100, 1, 20, 50, 10, 20))
#'
#' comat=df_to_comat(df,weight=TRUE)
#' @export
df_to_comat <- function(df, weight = FALSE){

  # Control
  if(!is.data.frame(df)){
    stop("df must be a data.frame with columns species and sites.")
  }

  if(dim(df)[2]!=2 & dim(df)[2]!=3){
    stop("df must be a two- or three-columns data.frame")
  }

  if(weight & dim(df)[2]==2){
    stop("df must be a three-columns data.frame if weight equal TRUE")
  }

  if(weight & dim(df)[2]==3 & class(df[,3])!="numeric" & class(df[,3])!="integer"){
    stop("The third column must be numeric")
  }

  # Rename columns
  colnames(df)[1:2]=c("Site","Species")
  if(weight){
    colnames(df)[3] <- "Weight"
  } else{ # without abundances, species are just present in the sites
    df$Weight <- 1
  }

  # Create contingency table with matrix indexing
  df$Site <- as.factor(df$Site)
  df$Species <- as.factor(df$Species)

  comat <- with(df, {
    out <- matrix(nrow = nlevels(Site), ncol = nlevels(Species),
                  dimnames = list(levels(Site), levels(Species)))
    out[cbind(Site, Species)] <- Weight
    out
  })
  # Replace NAs with 0s
  comat[is.na(comat)] <- 0

  # Check for empty rows and columns
  comat <- comat[rowSums(comat) > 0, ]
  comat <- comat[, colSums(comat) > 0]

  return(comat)
}
