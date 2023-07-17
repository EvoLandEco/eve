edd_summarize <- function(raw_data, method = NULL) {
  nrep <- length(raw_data$tes)
  statistics <- data.frame(id = 1:nrep)

  colless <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "Colless")
  j_one <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "J-One")
  tci <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "TCI")
  steps <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "Steps")
  branch <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "Branch")
  b1 <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "B1")
  b2 <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "B2")
  gamma <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "Gamma")
  mbl <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "MBL")
  pd <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "PD")
  mntd <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "MNTD")
  mpd <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "MPD")
  rogers <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "Rogers")
  sr <- lapply(raw_data$tes, calculate_tree_stats, method = method, metric = "SR")

  result <- cbind(Colless = unlist(colless), J_One = unlist(j_one), TCI = unlist(tci),
                  Steps = unlist(steps), Branch = unlist(branch), B1 = unlist(b1), B2 = unlist(b2),
                  Gamma = unlist(gamma), MBL = unlist(mbl), PD = unlist(pd),
                  MNTD = unlist(mntd), MPD = unlist(mpd), Rogers = unlist(rogers), SR = unlist(sr))

  statistics <- cbind(statistics, result)

  return(within(statistics, rm(id)))
}


#' @title Function to calculate tree statistics of search_tree results
#' @param data A list of trees
#' @param combo A list of parameter combinations
#' @return a data frame of tree statistics
#' @author Tianjian Qin
#' @export edd_stat_rtree
edd_stat_rtree <- function(data, combo) {
  data <- lapply(data, function (x) {
    x <- list(tes = x)
    return(x)
  })

  params <- combo_to_tibble(combo)
  meta <- params[rep(seq_len(nrow(params)), each = length(data[[1]]$tes)), ]
  statistics <- lapply(data, edd_summarize, method = "treestats")

  stats <- cbind(meta, dplyr::bind_rows(statistics))

  return(stats)
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


#' @title Function to calculate tree statistics of simulation results
#' @param raw_data Results of simulation
#' @param method Method to calculate tree statistics, default is "treestats"
#' @param strategy Parallel strategy, default is "sequential"
#' @param workers Number of workers, default is 1
#' @param verbose Whether to print progress, default is TRUE
#' @return a data frame of tree statistics
#' @author Tianjian Qin
#' @export edd_stat
edd_stat <- function(raw_data, method = "treestats", strategy = "sequential",
                     workers = 1, verbose = TRUE) {
  check_parallel_arguments(strategy, workers, verbose)
  check_raw_data(raw_data)

  statistics <- furrr::future_map(.x = raw_data,
                                  .f =  edd_summarize,
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


edd_stat_phylogenetic_evenness <- function(raw_data = NULL) {
  nrep <- length(raw_data$tes)
  statistics <- data.frame(id = 1:nrep)

  ERE <- mapply(calculate_phylogenetic_evenness, raw_data$tes, raw_data$las)
  result <- data.frame(ERE = unlist(ERE))
  statistics <- cbind(statistics, result)

  return(within(statistics, rm(id)))
}


edd_phylogenetic_evenness <- function(raw_data = NULL, strategy = "sequential",
                              workers = 1, verbose = TRUE) {
  check_parallel_arguments(strategy, workers, verbose)
  check_raw_data(raw_data)

  statistics <- furrr::future_map(.x = raw_data,
                                  .f =  edd_stat_phylogenetic_evenness)

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