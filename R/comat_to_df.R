#' Create a data.frame from a contingency table
#'
#' This function creates a two- or three-columns data.frame where
#' each row represents the presence of a species (second column)
#' in a site (first column) and an optional third column indicating the abundance
#' of the species in this site (if weight=TRUE) from a co-occurence matrix
#' (sites as rows and species as columns).
#'
#' @param comat a contingency matrix with sites as rows and species as columns
#' @param weight a boolean indicating if the value are presence (i.e. equal 1) or abundance data
#' @param remove_absent_species a boolean determining whether absent
#' species from the co-occurence matrix have to be removed from the output.
#' @export
#' @return A data.frame  where each row represents the presence of a species (second column)
#' in a site (first column) and an optional third column indicating the abundance
#' of the species in this site
#' @author
#' Pierre Denelle \email{pierre.denelle@gmail.com}
#' Maxime Lenormand \email{maxime.lenormand@inrae.fr}
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' comat=matrix(sample(1000,50),5,10)
#' rownames(comat)=paste0("Site",1:5)
#' colnames(comat)=paste0("Species",1:10)
#'
#' df=comat_to_df(comat,weight=TRUE)
#' @export
comat_to_df <- function(comat, weight = FALSE, remove_absent_species = TRUE){

  if(!is.matrix(comat)){
    stop("Contingency table should be a matrix with sites as rows and species as columns.")
  }

  if(!is.logical(remove_absent_species)){
    stop("remove_absent_species must be a boolean determining whether absent
    species from the contingency matrix have to be removed from the output.")
  }

  # Conversion as data.frame
  df <- reshape2::melt(codat)
  colnames(df) <- c("Site", "Species", "Weight")

  if(remove_absent_species == TRUE){
    df <- df[df$Weight > 0,]
  }

  if(!weight){
    df=df[,-3]
  }

  return(df)
}
