edd_summarize <- function(raw_data, stat = NULL, method = NULL) {
  nrep <- length(raw_data$tes)
  statistics <- data.frame(id = 1:nrep)

  if ("all" %in% stat | "balance" %in% stat) {
    sackin <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Sackin")
    colless <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Colless")
    blum <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Blum")
    balance <- cbind(sackin = unlist(sackin), colless = unlist(colless), blum = unlist(blum))

    statistics <- cbind(statistics, balance)
  }

  if ("all" %in% stat | "mbl" %in% stat) {
    mean_branch_length <- lapply(raw_data$tes, calculate_mean_branch_length, method = method)
    statistics <- cbind(statistics, mbl = unlist(mean_branch_length))
  }

  if ("all" %in% stat | "gamma" %in% stat) {
    gamma_statistics <- lapply(raw_data$tes, calculate_gamma_statistics, method = method)
    statistics <- cbind(statistics, gamma = unlist(gamma_statistics))
  }

  if ("all" %in% stat | "pd" %in% stat) {
    phylogenetic_diversity <- lapply(raw_data$tes, calculate_phylogenetic_diversity, method = method)
    statistics <- cbind(statistics, pd = unlist(phylogenetic_diversity))
  }

  if ("all" %in% stat | "mntd" %in% stat) {
    mean_nearest_neighbor_distance <- lapply(raw_data$tes, calculate_mean_nearest_neighbor_distance, method = method)
    statistics <- cbind(statistics, mntd = unlist(mean_nearest_neighbor_distance))
  }

  return(within(statistics, rm(id)))
}



edd_stat <- function(raw_data, stat = "all", method = "treestats", strategy = "sequential",
                     workers = 1, verbose = TRUE) {
  check_parallel_arguments(strategy, workers, verbose)
  check_raw_data(raw_data)

  statistics <- furrr::future_map(.x = raw_data,
                                  .f = edd_summarize,
                                  stat = stat,
                                  method = method)

  # Binding metadata
  nrep <- length(raw_data$`1`$tes)
  meta <- lapply(sim_data, function(x) {
    meta_data <- as.data.frame(extract_parameters(x))
    dplyr::slice(meta_data, rep(1, nrep))
  })
  meta <- dplyr::bind_rows(meta)
  statistics <- dplyr::bind_rows(statistics)
  statistics <- cbind(meta, statistics)

  return(statistics)
}
