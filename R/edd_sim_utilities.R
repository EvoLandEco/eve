#' @name L2ED_cpp
#' @title C++ version of L2ED
#' @description Function to convert a table with speciation and extinction events to
#' mean evolutionary distances between each species and the rest of the community
#' @param L Matrix of events as produced by edd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Simulation time
#' @return a named vector of mean evolutionary distinctivenesses
#' @author Tianjian Qin
#' @export L2ED_cpp
L2ED_cpp <- function(L, t) {
  dist_tips <-
    ape::cophenetic.phylo(treestats::l_to_phylo_ed(L, t, drop_extinct = TRUE))
  dist_means <- rowSums(dist_tips) / (dim(dist_tips)[1] - 1)
  dist_means_sorted <-
    dist_means[gtools::mixedorder(names(dist_means))]

  return(dist_means_sorted)
}



#' @name L2NND
#' @title Converting a table with speciation and extinction events to nearest neighbor distances
#' @description Function to convert a table with speciation and extinction events to
#' phylogenetic distances between each species and its nearest neighbor on the tree
#' @param L Matrix of events as produced by edd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Simulation time
#' @return a named vector of mean nearest neighbor distances
#' @author Tianjian Qin
#' @export L2NND
L2NND <- function(L, t) {
  dist_tips <-
    ape::cophenetic.phylo(treestats::l_to_phylo_ed(L, t, drop_extinct = TRUE))
  # find the second smallest value of each row of the distance matrix
  nearest_neighbor <- apply(dist_tips, MARGIN = 1, FUN = Rfast::nth, 2)
  nearest_neighbor_sorted <-
    nearest_neighbor[gtools::mixedorder(names(nearest_neighbor))]

  return(nearest_neighbor_sorted)
}



#' @name L2ED
#' @title Converting a table with speciation and extinction events to evolutionary
#' distances
#' @description Function to convert a table with speciation and extinction events to
#' mean evolutionary distances between each species and the rest of the community
#' @param L Matrix of events as produced by pdd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Simulation time
#' @return a named vector of mean evolutionary distinctivenesses
#' @author Tianjian Qin
#' @keywords models
#' @export L2ED
L2ED <- function(L, t) {
  dist_tips <-
    ape::cophenetic.phylo(L2phylo2(L, t, dropextinct = TRUE))
  dist_means <- rowSums(dist_tips) / (dim(dist_tips)[1] - 1)
  dist_means_sorted <-
    dist_means[gtools::mixedorder(names(dist_means))]

  return(dist_means_sorted)
}



#' @name L2Phi_cpp
#' @title C++ version of L2Phi
#' #' @description Function to convert a table with speciation and extinction events to a
#' phylogenetic diversity metric
#' @param L Matrix of events as produced by pdd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Sets whether the phylogeny should drop species that are
#' extinct at the present
#' @param metric Specifies which phylogenetic diversity metric should be used
#' @return a value of one of the phylogenetic diversity metrices
#' @author Tianjian Qin
#' @keywords models
#' @export L2Phi_cpp
L2Phi_cpp <- function(L, t, metric) {
  # metrics
  if (metric == "pd") {
    return(sum(treestats::l_to_phylo_ed(L, t, drop_extinct = T)$edge.length))
  } else if (metric == "mpd") {
    phy <- treestats::l_to_phylo_ed(L, t, drop_extinct = T)
    n <- length(phy$tip.label)
    dist <- ape::dist.nodes(phy)[1:n, 1:n]
    return(mean(dist[lower.tri(dist)]))
  }
}



#' @name L2Phi
#' @title Converting a table with speciation and extinction events to a phylogenetic
#' diversity metric
#' @description Function to convert a table with speciation and extinction events to a
#' phylogenetic diversity metric
#' @param L Matrix of events as produced by pdd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Sets whether the phylogeny should drop species that are
#' extinct at the present
#' @param metric Specifies which phylogenetic diversity metric should be used
#' @return a value of one of the phylogenetic diversity metrices
#' @author Tianjian Qin
#' @keywords models
#' @export L2Phi
L2Phi <- function(L, t, metric) {
  # metrics
  if (metric == "pd") {
    return(sum(L2phylo2(L, t, dropextinct = T)$edge.length))
  } else if (metric == "mpd") {
    phy <- L2phylo2(L, t, dropextinct = T)
    n <- length(phy$tip.label)
    dist <- ape::dist.nodes(phy)[1:n, 1:n]
    return(mean(dist[lower.tri(dist)]))
  }
}



#' Function to convert a table with speciation and extinction events to a
#' phylogeny
#'
#' Converting a table with speciation and extinction events to a phylogeny
#' (reversed time scale)
#'
#'
#' @param L Matrix of events as produced by dd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Simulation time
#' @param dropextinct Sets whether the phylogeny should drop species that are
#' extinct at the present
#' @return \item{ phy }{ A phylogeny of the phylo type }
#' @author Rampal S. Etienne; Tianjian Qin
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' # do not use this function, use L2phylo() instead
#'
#' @export L2phylo2
L2phylo2 <- function(L, t, dropextinct = FALSE)
  # makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
  L2 <- L[order(abs(L[, 3])), 1:4]
  age <- t
  L2[1, 1] <- -1
  sall <- which(L2[, 4] >= -1)
  tend <- (L2[, 4] == -1) * age + (L2[, 4] > -1) * L2[, 4]
  L2 <- L2[, -4]
  linlist <-
    cbind(data.frame(L2[sall,]), paste0("t", abs(L2[sall, 3])), tend)
  linlist[, 4] <- as.character(linlist[, 4])
  names(linlist) <- 1:5
  done <- 0
  while (done == 0) {
    j <- which.max(linlist[, 1])
    parent <- linlist[j, 2]
    parentj <- which(parent == linlist[, 3])
    parentinlist <- length(parentj)
    if (parentinlist == 1) {
      spec1 <-
        paste0(linlist[parentj, 4], ":", linlist[parentj, 5] - linlist[j, 1])
      spec2 <-
        paste0(linlist[j, 4], ":", linlist[j, 5] - linlist[j, 1])
      linlist[parentj, 4] <-
        paste0("(", spec1, ",", spec2, ")")
      linlist[parentj, 5] <- linlist[j, 1]
      linlist <- linlist[-j,]
    } else {
      linlist[j, 1:3] <- L2[which(L2[, 3] == parent), 1:3]
    }
    if (nrow(linlist) == 1) {
      done <- 1
    }
  }
  linlist[4] <- paste0(linlist[4], ":", linlist[5], ";")
  phy <- ape::read.tree(text = linlist[1, 4])
  tree <- ape::as.phylo(phy)

  if (dropextinct == FALSE) {
    return(tree)
  } else {
    absent <- abs(L[which(L[, 4] != -1), 3])
    tree <- L2phylo2(L, t, dropextinct = F)
    pruned_tree <-
      ape::drop.tip(tree, paste0("t", absent))

    return(pruned_tree)
  }
}




#' @name phylo2mpd
#' @title Converting a phylogenetic tree of phylo object to mean pairwise distance
#' among all tips
#' @description Function to convert a phylogenetic tree with extant species to mean pairwise
#' distance
#' @param phy a phylo object that only contains extant species
#' extinct at the present
#' @return a numeric value of mean pairwise distance
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export phylo2mpd
phylo2mpd <- function(phy) {
  n <- length(phy$tip.label)
  dist <- ape::dist.nodes(phy)[1:n, 1:n]
  return(mean(dist[lower.tri(dist)]))
}



#' @name phylo2pd
#' @title Converting a phylogenetic tree of phylo object to total lengths of all
#' branches in the tree
#' @description Function to convert a phylogenetic tree with extant species to
#' phylogenetic diversity
#' @param phy a phylo object that only contains extant species
#' extinct at the present
#' @return a numeric value of phylogenetic diversity
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export phylo2pd
phylo2pd <- function(phy) {
  return(sum(phy$edge.length))
}



#' Takes samples in the usual manner
#'
#' The standard sample function in R samples from n numbers when x = n. This is
#' unwanted behavior when the size of the vector to sample from changes
#' dynamically. This is corrected in sample2
#'
#'
#' @param x A vector of one or more elements
#' @param size A non-negative integer giving the number of items to choose.
#' @param replace Should sampling be with replacement?
#' @param prob A vector of probability weights for obtaining the elements of
#' the vector being sampled.
#' @return \item{sam}{A vector of length \code{size} that is sampled from
#' \code{x}. }
#' @author Rampal S. Etienne
#' @keywords models
#' @examples
#'
#' sample(x = 10,size = 5,replace = TRUE)
#' sample2(x = 10,size = 5,replace = TRUE)
#'
#' @export sample2
sample2 <-  function(x,size,replace = FALSE,prob = NULL)
{
  if(length(x) == 1)
  {
    x <- c(x,x)
    prob <-  c(prob,prob)
    if(is.null(size))
    {
      size <-  1
    }
    if(replace == FALSE & size > 1)
    {
      stop('It is not possible to sample without replacement multiple times from a single item.')
    }
  }
  sam <-  sample(x,size,replace,prob)
  return(sam)
}
