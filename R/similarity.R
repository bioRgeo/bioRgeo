
similarity <- function(dat, similarity = "simpson", input = "matrix",
                            site = NULL, sp = NULL, ab = NULL, binary = TRUE){
  require(Rcpp)
  require(SMUT)

  if(!(input %in% c("matrix", "data frame"))){
    stop("dat must be either a contingency matrix with sites as rows and
    species as columns or a long format data frame with each row being the
    presence of a species in a given site.")
  }

  if(input == "matrix"){
    if(!is.matrix(dat)){
      stop("dat should be a matrix with sites as rows and species as columns.")
    }
  } else if(input == "data frame"){
    if(!is.data.frame(dat)){
      stop("dat should be a long format data frame with each row being the
    presence of a species in a given site.")
    }

    # Conversion as a contingency table with contingency function
    dat <- contingency(dat, site, sp, ab = NULL, binary = TRUE)
  }

  if(!(similarity %in% c("simpson", "jaccard", "sorensen", "whittaker",
                         "bray", "euclidean"))){
    stop("Similarity metric chosen is not available.
     Please chose among the followings:
         simpson, jaccard, sorensen, whittaker, bray or euclidean")
  }

  # Convert as matrix to remove names
  dat <- as.matrix(dat)

  if(similarity == "bray"){
    require(ecodist)
    require(reshape2)
    bray <- bcdist(dat, rmzero = FALSE)
    # Conversion as matrix and removal of upper part and diagonal
    bray <- as.matrix(bray)
    bray[upper.tri(bray)] <- NA
    diag(bray) <- NA

    # Convert distance matrix into data.frame
    abc <- reshape2::melt(bray, varnames = c("id1", "id2"))
    # Remove NAs
    abc <- abc[complete.cases(abc), ]
    colnames(abc) <- c("id1", "id2", "bray")
    # Create columns of site names and put id in the id columns
    if(!is.null(rownames(dat))){
      abc$id1_name <- abc$id1
      abc$id2_name <- abc$id2
      abc$id1 <- match(abc$id1, rownames(dat))
      abc$id2 <- match(abc$id2, rownames(dat))
    }
    # Similarity instead of dissimilarity
    abc$bray <- 1 - abc$bray

    # Remove zeros
    abc <- abc[which(abc$bray > 0), ]

  } else if(similarity == "euclidean"){
    require(reshape2)
    euc_dist <- function(m) {
      mtm <- Matrix::tcrossprod(m)
      sq <- rowSums(m*m)
      sqrt(outer(sq, sq, "+") - 2*mtm)
    }

    abc <- euc_dist(dat)
    # removal of upper part and diagonal
    abc[upper.tri(abc)] <- NA
    diag(abc) <- NA

    # Conversion to dataframe
    abc <- reshape2::melt(abc, varnames = c("id1", "id2"))
    # Remove NAs
    abc <- abc[complete.cases(abc), ]
    colnames(abc) <- c("id1", "id2", "euclid")

  } else{
    # Get the number of species shared by two sites
    # eigenMapMatMult to compute matrix product; equivalent to dat %*% t(dat)
    a <- eigenMapMatMult(dat, t(dat))

    abc <- which(a > 0, arr.ind = TRUE)
    abc <- cbind(abc, a[a > 0])
    # Remove upper part of the matrix
    abc <- abc[abc[, 1] < abc[, 2], ]

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
    if(similarity == "simpson"){
      abc$simpson <- 1 - pmin(abc$b, abc$c)/(abc$a + pmin(abc$b, abc$c))
      abc <- abc[, c("id1", "id2", "simpson")]
    } else if(similarity == "sorensen"){
      abc$sorensen <- 1 - (abc$b + abc$c)/(2*abc$a + abc$b + abc$c)
      abc <- abc[, c("id1", "id2", "sorensen")]
    } else if(similarity == "jaccard"){
      abc$jaccard <- 1 - abc$a/(abc$a + abc$b + abc$c)
      abc <- abc[, c("id1", "id2", "jaccard")]
    } else if(similarity == "whittaker"){
      abc$whittaker <- 1 - (abc$a + abc$b + abc$c)/((2*abc$a + abc$b + abc$c)/2)
      abc <- abc[, c("id1", "id2", "whittaker")]
    }
  }
  # If contingency matrix has rownames, reassign them to abc data.frame
  if(!is.null(rownames(dat))){
    abc$id1 <- rownames(dat)[abc$id1]
    abc$id2 <- rownames(dat)[abc$id2]
  }

  return(abc)
}
