# Method from is.extinct in package Geiger (Pennell, M.W., J.M. Eastman, G.J. Slater, J.W. Brown, J.C. Uyeda,
# R.G. FitzJohn, M.E. Alfaro, and L.J. Harmon. 2014. geiger v2.0: an expanded suite of methods for fitting
# macroevolutionary models to phylogenetic trees. Bioinformatics 30:2216-2218.) for the following function,
# with improvements and enhancement.
# By default returns a vector of edge numbers relative to the extinct tips, but can also return a vector of
# extinct tip labels.
find_extinct <- function(phy, tolerance = NULL, which = "edge_number") {
  if (!"phylo" %in% class(phy)) {
    stop("\"phy\" is not of class \"phylo\".")
  }
  if (is.null(phy$edge.length)) {
    stop("\"phy\" does not have branch lengths.")
  }
  if (is.null(tolerance)) {
    tolerance <- min(phy$edge.length) / 100
  }
  phy <- reorder(phy)
  tips <- numeric(ape::Ntip(phy) + phy$Nnode)
  for (i in seq_along(phy$edge[, 1])) {
    tips[phy$edge[i, 2]] <- tips[phy$edge[i, 1]] + phy$edge.length[i]
  }
  extinct_tips <- max(tips[1:ape::Ntip(phy)]) - tips[1:ape::Ntip(phy)] > tolerance

  if (any(extinct_tips)) {
    if (which == "tip_label") {
      return(phy$tip.label[which(extinct_tips)])
    } else if (which == "edge_number") {
      edge_number <- sapply(phy$tip.label[which(extinct_tips)], USE.NAMES = FALSE, function(x) {
        ape::which.edge(phy, x)
      })
      return(edge_number)
    } else {
      stop("Invalid value for \"which\"")
    }
  } else {
    return(NULL)
  }
}



mark_extinct_tips <- function(phy) {
  # Mark extinct tips with different linetype
  tree_data <- as(phy, "phylo4")
  extinct_tips <- find_extinct(phy, which = "tip_label")
  extinct_data <- data.frame(linetype = rep(1, length(phy$tip.label)))
  rownames(extinct_data) <- phy$tip.label
  extinct_data[extinct_tips, ] <- 2
  extinct_data$linetype <- as.factor(extinct_data$linetype)
  tree_data <- phylobase::phylo4d(tree_data, extinct_data)
  internal_node_data <- data.frame(linetype = rep(1, phylobase::nNodes(tree_data)),
                                   row.names = phylobase::nodeId(tree_data, "internal"))
  phylobase::nodeData(tree_data) <- internal_node_data

  return(tree_data)
}