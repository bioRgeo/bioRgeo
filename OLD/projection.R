
simil <- function(dat, metric = "simpson", input = "matrix",
                  output = "data.frame",
                  site = NULL, sp = NULL, ab = NULL, weight = FALSE){
  ## Controls ----
  require(Rcpp)
  require(SMUT)

  if(!(input %in% c("matrix", "data.frame"))){
    stop("dat must be either a contingency matrix with sites as rows and
    species as columns or a long format data.frame with each row being the
    presence of a species in a given site.")
  }

  if(!(output %in% c("matrix", "data.frame", "dist"))){
    stop("output is a character string indicating the format of the output.
         Either 'matrix', 'data.frame' or 'dist'.")
  }

  if(input == "matrix"){
    if(!is.matrix(dat)){
      stop("dat should be a matrix with sites as rows and species as columns.")
    }
  } else if(input == "data.frame"){
    if(!is.data.frame(dat)){
      stop("dat should be a long format data.frame with each row being the
    presence (or abundance) of a species in a given site.")
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
      warning("Without column abundances, contingency table will only get
      binary values.")
    }

    if(!is.logical(weight)){
      stop("weight must be a boolean.")
    }

    if(!is.null(ab) & metric != "bray"){
      warning("Only Bray-Curtis similarity metric can handle abundances of
              species within sites. Other metrics deal with presence/absence
              data only.")
    }

    # Conversion as a contingency table with contingency function
    dat <- contingency(dat, site, sp, ab = ab)
  }

  if(!(metric %in% c("simpson", "jaccard", "sorensen", "whittaker", "bray",
                     "euclidean"))){
    stop("Similarity metric chosen is not available.
     Please chose among the followings:
         simpson, jaccard, sorensen, whittaker, bray or euclidean")
  }

  ## Function ----
  # Convert as matrix to remove names
  dat <- as.matrix(dat)

  if(metric == "bray"){
    require(ecodist)
    require(reshape2)
    bray <- bcdist(dat, rmzero = FALSE)
    # Conversion as matrix and removal of upper part and diagonal
    bray <- as.matrix(bray)

    if(output == "dist"){
      # Similarity instead of dissimilarity
      bray <- 1 - bray
      # Conversion to a dist object
      abc <- as.dist(bray)
    }
    else{
      # Removal of upper part and diagonal of the matrix
      bray[upper.tri(bray)] <- NA
      diag(bray) <- NA

      # Similarity instead of dissimilarity
      bray <- 1 - bray

      if(output == "data.frame"){
        # Convert distance matrix into data.frame
        abc <- reshape2::melt(bray, varnames = c("id1", "id2"))
        # Remove NAs
        abc <- abc[complete.cases(abc), ]
        colnames(abc) <- c("id1", "id2", "bray")

        # Remove zeros
        # abc <- abc[which(abc$bray > 0), ]

        # If contingency matrix has rownames, reassign them to abc data.frame
        if(output == "data.frame" & !is.null(rownames(dat))){
          abc$id1 <- rownames(dat)[abc$id1]
          abc$id2 <- rownames(dat)[abc$id2]
        }
      }
    }
  } else if(metric == "euclidean"){
    require(reshape2)
    euc_dist <- function(m) {
      mtm <- Matrix::tcrossprod(m)
      sq <- rowSums(m*m)
      sqrt(outer(sq, sq, "+") - 2*mtm)
    }

    abc <- euc_dist(dat)

    if(output == "dist"){
      abc <- as.dist(abc)
    }

    if(output == "data.frame"){
      # Removal of upper part and diagonal
      abc[upper.tri(abc)] <- NA
      diag(abc) <- NA

      # Conversion to dataframe
      abc <- reshape2::melt(abc, varnames = c("id1", "id2"))
      # Remove NAs
      abc <- abc[complete.cases(abc), ]
      colnames(abc) <- c("id1", "id2", "euclid")

      # If contingency matrix has rownames, reassign them to abc data.frame
      if(output == "data.frame" & !is.null(rownames(dat))){
        abc$id1 <- rownames(dat)[abc$id1]
        abc$id2 <- rownames(dat)[abc$id2]
      }
    }

  } else{
    # Get the number of species shared by two sites
    # eigenMapMatMult to compute matrix product; equivalent to dat %*% t(dat)
    a <- eigenMapMatMult(dat, t(dat))

    abc <- which(a > 0, arr.ind = TRUE)
    abc <- cbind(abc, a[a > 0])

    # Sites entirely dissimilar (no species shared)
    abc_0 <- which(a == 0, arr.ind = TRUE)
    abc_0 <- cbind(abc_0, a[a == 0])

    abc <- rbind(abc, abc_0)
    rm(abc_0)

    if(output == "data.frame"){
      # Remove upper part of the matrix (not if outputs = 'matrix' or 'dist' to
      # have symmetrical matrices)
      abc <- abc[abc[, 1] < abc[, 2], ]
    }

    a <- diag(a) # number of species per sites
    b <- a[abc[, 1]] # number of species in first column
    c <- a[abc[, 2]] # number of species in second column

    # Bind informations
    abc <- cbind(abc, b, c)
    # b equals to species only in first column
    abc[, 4] <- abc[, 4] - abc[, 3]
    # c equals to species only in second column
    abc[, 5] <- abc[, 5] - abc[, 3]

    colnames(abc) <- c("id1", "id2", "a", "b", "c")

    # Conversion as data.frame
    abc <- as.data.frame(abc)

    # Similarity metric
    if(metric == "simpson"){
      abc$simpson <- 1 - pmin(abc$b, abc$c)/(abc$a + pmin(abc$b, abc$c))
      abc <- abc[, c("id1", "id2", "simpson")]
    } else if(metric == "sorensen"){
      abc$sorensen <- 1 - (abc$b + abc$c)/(2*abc$a + abc$b + abc$c)
      abc <- abc[, c("id1", "id2", "sorensen")]
    } else if(metric == "jaccard"){
      abc$jaccard <- abc$a/(abc$a + abc$b + abc$c)
      abc <- abc[, c("id1", "id2", "jaccard")]
    } else if(metric == "whittaker"){
      abc$whittaker <- 1 - (abc$a + abc$b + abc$c)/((2*abc$a + abc$b + abc$c)/2)
      abc <- abc[, c("id1", "id2", "whittaker")]
    }

    # If contingency matrix has rownames, reassign them to abc data.frame
    if(!is.null(rownames(dat))){
      abc$id1 <- rownames(dat)[abc$id1]
      abc$id2 <- rownames(dat)[abc$id2]
    }

    if(output %in% c("matrix", "dist")){
      require(reshape2)
      abc <- reshape2::dcast(abc, id1 ~ id2,
                             value.var = colnames(abc)[3])
      rownames(abc) <- abc$id1
      abc <- abc[, !(colnames(abc) == "id1")]
      if(output == "dist"){
        abc <- as.dist(as.matrix(abc))
      }
    }
  }
  return(abc)
}
