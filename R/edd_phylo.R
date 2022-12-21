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



sample_end_state <- function(l_table, params, which = stop("Please specify which end state to sample")) {
  metric <- params$metric
  model <- params$model
  age <- params$age
  num <- nrow(l_table[l_table[, 4] == -1,])
  linlist <- l_table[l_table[, 4] == -1,][, 3]
  if (model == "dsce2") {
    pars <- c(num, unlist(params$pars)[1],
              unlist(params$pars)[2],
              unlist(params$pars)[3],
              unlist(params$pars)[4])
  }
  if (model == "dsde2") {
    pars <- c(num, unlist(params$pars)[1],
              unlist(params$pars)[2],
              unlist(params$pars)[3],
              unlist(params$pars)[4],
              unlist(params$pars)[5],
              unlist(params$pars)[6])
  }

  if (metric == "pd") {
    ed <- rep(as.vector(L2Phi_cpp(l_table, age, "pd")), num)
  }
  if (metric == "ed") {
    ed <- L2ED_cpp(l_table, age)
  }
  if (metric == "nnd") {
    ed <- L2NND(l_table, age)
  }

  rates <- edd_update_lamu(ed, ed, pars, model)
  las <- purrr::set_names(rates$newlas, paste0('t', abs(linlist)))
  las <- c(time = age, las)
  mus <- purrr::set_names(rates$newmus, paste0('t', abs(linlist)))
  mus <- c(time = age, mus)
  eds <- purrr::set_names(ed, paste0('t', abs(linlist)))
  eds <- c(time = age, eds)

  if (which == "las") {
    return(las)
  } else if (which == "mus") {
    return(mus)
  } else if (which == "eds") {
    return(eds)
  }

  return(list(las = las, mus = mus, eds = eds))
}


#' @name stat_histree
#' @title Mapping historical states to a phylo tree
#' @description Function to map historical states to a phylo tree and output a data frame of segments of the tree
#' @param phy A phylo object, the tree to map historical states to
#' @param history A data frame containing time steps and the historical states of each lineage. The first row of the data
#' frame should be time, which has each time step as a row. The other rows should be named t1, t2, t3, etc, and the names
#' should be in the same order as they speciated. The rows should contain the values of its historical states at each time
#' step. At each time step where a lineage is not present, the value should be NA.
#' @return A data frame of the segments of the tree
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export stat_histree
stat_histree <- function(phy, history){
  # get the skeleton of the phylo tree, get all the key coordinates
  xs <- ape::node.depth.edgelength(phy)
  ys <- ape::node.height(phy)
  nodes <- data.frame(x = xs, y = ys)

  # find the root node and its coordinate
  root_node <- nodes[which(nodes$x == 0),]
  root_node$id <- ape::Ntip(phy) + 1

  segments <- data.frame()
  # an iterator function to look up and draw the state history of two daughter edges connected to a node
  draw_daughter_lineages(parent_node =  root_node,
                         phy = phy,
                         history= history,
                         nodes = nodes,
                         tolerance = 1e-8)

  return(segments)
}

#' @name draw_daughter_lineages
#' @title Draw the state history of two daughter edges connected to a node
#' @description An iterator function to look up and draw the state history of two daughter edges connected to a node
draw_daughter_lineages <- function(envir = parent.frame(), parent_node, phy, history, nodes, tolerance = 1e-9) {
  tolerance_r <- 1 / tolerance
  parent_age <- parent_node$x
  daughters <- phy$edge[phy$edge[, 1] == parent_node$id,]
  daughters_nodes <- daughters[, 2]
  daughters_lengths <- phy$edge.length[which(phy$edge[, 1] == parent_node$id)]

  # find the representative descendant of each daughter lineage, which is the earliest speciated descendant,
  # note that this function requires the lineages to be label like t1, t2, ... tn, and in the order of the
  # time of speciation, that's why gtools::mixedsort() is required here, after all, the earliest speciated
  # descendant has the full history to be able to map onto the shared ancestral edges
  # I mis-spelled descendant here and after
  decendants1 <- geiger::tips(phy, daughters_nodes[1])


  # find matching history given the time frame (age of parent node plus edge length)
  if (length(decendants1) > 1) {
    rep_decendants1 <- gtools::mixedsort(decendants1)[1]
    history1 <- history %>%
      dplyr::select(time, rep_decendants1) %>%
      dplyr::filter(time >= (floor(parent_age * tolerance_r) / tolerance_r) &
                      time <= (ceiling((parent_age + daughters_lengths[1]) * tolerance_r) / tolerance_r))
    daughter1_coord <- nodes[dplyr::near(nodes$x, parent_age + daughters_lengths[1], tolerance),]
  } else {
    rep_decendants1 <- decendants1
    history1 <- history %>%
      dplyr::select(time, rep_decendants1) %>%
      dplyr::filter(time >= (floor(parent_age * tolerance_r) / tolerance_r) &
                      time <= (ceiling((parent_age + daughters_lengths[1]) * tolerance_r) / tolerance_r))
    daughter1_coord <- nodes[which(phy$tip.label == rep_decendants1),]
  }

  # generate coordinates to arrange the horizontal segments
  segments1h <- data.frame()
  for (i in 1:(nrow(history1) - 1)) {
    x <- history1$time[i]
    xend <- history1$time[i + 1]
    y <- daughter1_coord$y
    yend <- daughter1_coord$y
    state <- (history1[rep_decendants1][i,] + history1[rep_decendants1][i + 1,]) / 2
    segments1hnew <- data.frame(x = x, y = y, xend = xend, yend = yend, state = state)
    segments1h <- rbind(segments1h, segments1hnew)
  }

  # generate coordinates for the two vertical segments (connected to the parent node)
  segments1v <- data.frame(x = parent_age,
                           y = parent_node$y,
                           xend = parent_age,
                           yend = daughter1_coord$y,
                           state = history1[rep_decendants1][1,])

  # repeat the same process for another descendant edge from the parent node
  decendants2 <- geiger::tips(phy, daughters_nodes[2])
  if (length(decendants2) > 1) {
    rep_decendants2 <- gtools::mixedsort(decendants2)[1]
    history2 <- history %>%
      dplyr::select(time, rep_decendants2) %>%
      dplyr::filter(time >= (floor(parent_age * tolerance_r) / tolerance_r) &
                      time <= (ceiling((parent_age + daughters_lengths[2]) * tolerance_r) / tolerance_r))
    daughter2_coord <- nodes[dplyr::near(nodes$x, parent_age + daughters_lengths[2], tolerance),]
  } else {
    rep_decendants2 <- decendants2
    history2 <- history %>%
      dplyr::select(time, rep_decendants2) %>%
      dplyr::filter(time >= (floor(parent_age * tolerance_r) / tolerance_r) &
                      time <= (ceiling((parent_age + daughters_lengths[2]) * tolerance_r) / tolerance_r))
    daughter2_coord <- nodes[which(phy$tip.label == rep_decendants2),]
  }

  segments2h <- data.frame()
  for (i in 1:(nrow(history2) - 1)) {
    x <- history2$time[i]
    xend <- history2$time[i + 1]
    y <- daughter2_coord$y
    yend <- daughter2_coord$y
    state <- (history2[rep_decendants2][i,] + history2[rep_decendants2][i + 1,]) / 2
    segments2hnew <- data.frame(x = x, y = y, xend = xend, yend = yend, state = state)
    segments2h <- rbind(segments2h, segments2hnew)
  }

  segments2v <- data.frame(x = parent_age,
                           y = parent_node$y,
                           xend = parent_age,
                           yend = daughter2_coord$y,
                           state = history2[rep_decendants2][1,])

  # amend segments data frame in the outer scope
  envir$segments <- rbind(envir$segments, segments1h, segments1v, segments2h, segments2v)

  # iteration starts here, apply the function itself to the end nodes of the two descendant edges
  if (length(decendants1) > 1) {
    daughter1 <- data.frame(x = daughter1_coord$x,
                            y = daughter1_coord$y,
                            id = daughters_nodes[1])
    if (nrow(daughter1) > 1) {
      stop("Error: daughter1 has more than one row")
    }
    draw_daughter_lineages(envir = envir,
                           parent_node = daughter1,
                           phy = phy,
                           history = history,
                           nodes = nodes)
  }

  if (length(decendants2) > 1) {
    daughter2 <- data.frame(x = daughter2_coord$x,
                            y = daughter2_coord$y,
                            id = daughters_nodes[2])
    if (nrow(daughter2) > 1) {
      stop("Error: daughter2 has more than one row")
    }
    draw_daughter_lineages(envir = envir,
                           parent_node = daughter2,
                           phy = phy,
                           history = history,
                           nodes = nodes)
  }
}


