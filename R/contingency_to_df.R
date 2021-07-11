#' Create a data.frame from a contingency table
#'
#' This function creates a two- or three-columns data.frame where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction (if weight=TRUE)
#' from a contingency table (sites as rows and species as columns for example).
#'
#' @param comat a contingency table
#' @param weight a boolean indicating if the value are weights
#' @param remove_absent_objects a boolean determining whether absent
#' objects from the contingency table have to be removed from the output
#' @export
#' @return A data.frame where each row represents the interaction between
#' two objects and an optional third column indicating the weight of the interaction
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
contingency_to_df <- function(comat, weight = FALSE, remove_absent_objects = TRUE){

  # Control
  if(!is.matrix(comat)){
    stop("Contingency table should be a matrix")
  }

  if(!is.logical(remove_absent_objects)){
    stop("remove_absent_objects must be a boolean")
  }

  # Conversion as data.frame
  df <- reshape2::melt(comat)
  colnames(df) <- c("Object 1", "Object 2", "Weight")

  if(remove_absent_objects == TRUE){
    df <- df[df$Weight > 0,]
  }

  if(!weight){
    df=df[,-3]
  }

  return(df)
}
