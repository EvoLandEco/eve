edd_summarize <- function(raw_data, method = NULL) {
  nrep <- length(raw_data$tes)
  statistics <- data.frame(id = 1:nrep)

  colless <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Colless")
  j_one <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "J-One")
  tci <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "TCI")
  steps <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Steps")
  branch <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Branch")
  b1 <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "B1")
  b2 <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "B2")
  gamma <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "Gamma")
  mbl <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "MBL")
  pd <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "PD")
  mntd <- lapply(raw_data$tes, calculate_tree_balance, method = method, metric = "MNTD")

  result <- cbind(Colless = unlist(colless), J_One = unlist(j_one), TCI = unlist(tci),
                   Steps = unlist(steps), Branch = unlist(branch), B1 = unlist(b1), B2 = unlist(b2),
                   Gamma = unlist(gamma), MBL = unlist(mbl), PD = unlist(pd), MNTD = unlist(mntd))

  statistics <- cbind(statistics, result)

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



edd_stat <- function(raw_data, method = "treestats", strategy = "sequential",
                     workers = 1, verbose = TRUE) {
  check_parallel_arguments(strategy, workers, verbose)
  check_raw_data(raw_data)

  edd_summarize_cached <- R.cache::addMemoization(edd_summarize)

  statistics <- furrr::future_map(.x = raw_data,
                                  .f =  edd_summarize_cached,
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