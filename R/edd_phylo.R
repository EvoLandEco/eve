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



# Reconstruct the growth of the phylo tree from the present L-table, with simulation time provided.
# Will output Blum, Colless, Sackin, PD, MNTD and MBL statistics for each time step.
reconstruct_temporal_dynamics <- function(l_table = NULL, age = NULL) {
  if (is.null(l_table)) {
      stop("No L-table provided")
  }
  if (is.null(age)) {
      stop("No age provided")
  }
  time_speciation <- unique(l_table[, 1])
  time_extinction <- unique(l_table[l_table[, 4] > 0, 4])
  event_table <- data.frame(Age = c(time_speciation, time_extinction),
                            Event = c(rep("speciation", length(time_speciation)),
                                      rep("extinction", length(time_extinction))))
  event_table <- event_table[order(event_table[, 1]),]
  temporal_dynamic <- data.frame(Age = 0,
                                 Event = "speciation",
                                 Blum = NA,
                                 Colless = NA,
                                 Sackin = NA,
                                 PD = 0,
                                 MNTD = 0,
                                 MBL = 0)

  # might be safe given often used rates
  moment <- 0.00001

  for (i in seq_along(event_table[, 1])) {
    if (i > 1) {
      if (event_table[i, 2] == "speciation") {
        t <- event_table[i, 1]
        l_temp <- l_table[l_table[, 1] < t,]
        l_temp[l_temp[, 4] > t, 4] <- -1
        temporal_dynamic <- rbind(temporal_dynamic,
                                  sample_tree_metrics(l_temp, t, event_table[i, 2], drop_extinct = TRUE))
      } else {
        # firstly sample the stats at the moment before extinction
        t <- event_table[i, 1] - moment
        l_temp <- l_table[l_table[, 1] < t,]
        l_temp[l_temp[, 4] > t, 4] <- -1
        temporal_dynamic <- rbind(temporal_dynamic,
                                  sample_tree_metrics(l_temp, t, event_table[i, 2], drop_extinct = TRUE))
        # then sample again the stats at the moment after extinction
        t <- event_table[i, 1] + moment
        l_temp <- l_table[l_table[, 1] < t,]
        l_temp[l_temp[, 4] > t, 4] <- -1
        temporal_dynamic <- rbind(temporal_dynamic,
                                  sample_tree_metrics(l_temp, t, event_table[i, 2], drop_extinct = TRUE))
      }
    }
  }

  temporal_dynamic <- rbind(temporal_dynamic,
                            sample_tree_metrics(l_table, age, "present", drop_extinct = TRUE))

  return(temporal_dynamic)
}



sample_tree_metrics <- function(l_table, age, event, drop_extinct) {
  phy <- treestats::l_to_phylo_ed(l_table, age, drop_extinct = drop_extinct)
  Blum <- calculate_tree_balance(phy, metric = "Blum")
  Colless <- calculate_tree_balance(phy, metric = "Colless")
  Sackin <- calculate_tree_balance(phy, metric = "Sackin")
  PD <- calculate_phylogenetic_diversity(phy)
  MNTD <- calculate_mean_nearest_neighbor_distance(phy)
  MBL <- calculate_mean_branch_length(phy)

  return(data.frame(Age = age,
                    Event = event,
                    Blum = Blum,
                    Colless = Colless,
                    Sackin = Sackin,
                    PD = PD,
                    MNTD = MNTD,
                    MBL = MBL))
}



