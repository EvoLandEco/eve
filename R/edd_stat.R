edd_summarize <- function(raw_data, stat = NULL, method = NULL) {
  nrep <- length(raw_data$tes)
  statistics <- data.frame(id = 1:nrep)

  if ("all" %in% stat | "balance" %in% stat) {
    sackin <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Sackin")
    colless <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Colless")
    blum <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Blum")
    balance <- cbind(Sackin = unlist(sackin), Colless = unlist(colless), Blum = unlist(blum))

    statistics <- cbind(statistics, balance)
  }

  if ("all" %in% stat | "mbl" %in% stat) {
    mean_branch_length <- lapply(raw_data$tes, calculate_mean_branch_length, method = method)
    statistics <- cbind(statistics, MBL = unlist(mean_branch_length))
  }

  if ("all" %in% stat | "gamma" %in% stat) {
    gamma_statistics <- lapply(raw_data$tes, calculate_gamma_statistics, method = method)
    statistics <- cbind(statistics, Gamma = unlist(gamma_statistics))
  }

  if ("all" %in% stat | "pd" %in% stat) {
    phylogenetic_diversity <- lapply(raw_data$tes, calculate_phylogenetic_diversity, method = method)
    statistics <- cbind(statistics, PD = unlist(phylogenetic_diversity))
  }

  if ("all" %in% stat | "mntd" %in% stat) {
    mean_nearest_neighbor_distance <- lapply(raw_data$tes, calculate_mean_nearest_neighbor_distance, method = method)
    statistics <- cbind(statistics, MNTD = unlist(mean_nearest_neighbor_distance))
  }

  if ("all" %in% stat | "sr" %in% stat) {
    tree_sizes <- get_tree_sizes(raw_data$tes)
    statistics <- cbind(statistics, SR = tree_sizes)
  }

  return(within(statistics, rm(id)))
}



edd_summarize_temporal_dynamics <- function(raw_data, rep_id = 1, strategy = "sequential",
                                            workers = 1, verbose = TRUE) {
  check_parallel_arguments(strategy, workers, verbose)
  check_raw_data(raw_data)

  temporal_dynamics <- furrr::future_map(.x = raw_data,
                                  .f = function (x) {
                                    stats <- reconstruct_temporal_dynamics(x$l_tables[[rep_id]], x$all_pars$age)
                                    meta <- as.data.frame(extract_parameters(x))
                                    return(cbind(meta, stats))
                                  })

  temporal_dynamics <- dplyr::bind_rows(temporal_dynamics)

  return(temporal_dynamics)
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
  meta <- lapply(raw_data, function(x) {
    meta_data <- as.data.frame(extract_parameters(x))
    dplyr::slice(meta_data, rep(1, nrep))
  })
  meta <- dplyr::bind_rows(meta)
  statistics <- dplyr::bind_rows(statistics)
  statistics <- cbind(meta, statistics)

  return(statistics)
}