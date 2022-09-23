check_parallel_arguments <- function(strategy = NULL,
                                     workers = NULL,
                                     verbose = TRUE) {
  if (workers == 1 |
      identical(strategy, "sequential")) {
    if (verbose == TRUE) {
      message("Running sequential task")
    }
    strategy <- eval(parse(text = paste0("future::", strategy)))
    future::plan(strategy)
  } else if (!(workers %% 1 == 0)) {
    stop("number of workers should be an integer")
  } else {
    if (strategy %in% c("multisession", "multicore", "cluster")) {
      if (verbose == TRUE) {
        message(paste0(
          "Running ",
          strategy,
          " task with ",
          workers,
          " workers"
        ))
      }
      strategy <- eval(parse(text = paste0("future::", strategy)))
      future::plan(strategy, workers = workers)
    } else {
      stop("incorrect parallel computing strategy")
    }
  }
}
