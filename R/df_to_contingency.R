#' Create a contingency table from a data.frame
#'
#' This function creates a contingency table from a two- or three-columns data.frame where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction (if weight=TRUE).
#'
#' @param df a two- or three-columns data.frame where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction
#' @param weight a boolean indicating if the abundance should be considered
#' @export
#' @return A contingency table with the first objects (first column of df) as rows and
#' the second objects (second column of df) as columns
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
df_to_contingency <- function(df, weight = FALSE){

  # Control
  if(!is.data.frame(df)){
    stop("df must be a two- or three-columns data.frame")
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
  colnames(df)[1:2]=c("Object1","Object2")

  # Manage weight
  if(weight){
    colnames(df)[3] <- "Weight"
  } else{
    df$Weight <- 1
  }

  # Create contingency table with matrix indexing
  df$Site <- as.factor(df$Object1)
  df$Species <- as.factor(df$Object2)

  comat <- with(df, {
    out <- matrix(nrow = nlevels(Object1), ncol = nlevels(Object2),
                  dimnames = list(levels(Object1), levels(Object2)))
    out[cbind(Object1, Object2)] <- Weight
    out
  })
  # Replace NAs with 0s
  comat[is.na(comat)] <- 0

  # Check for empty rows and columns
  comat <- comat[rowSums(comat) > 0, ]
  comat <- comat[, colSums(comat) > 0]

  return(comat)
}
