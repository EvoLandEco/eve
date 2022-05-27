check_parallel_arguments <- function(workers = NULL, strategy = NULL) {
  if (workers == 1 |
      identical(strategy, "sequential")) {
    message("Running sequential data loading")
    strategy <- eval(parse(text = paste0("future::", strategy)))
    future::plan(strategy)
  } else if (!(workers %% 1 == 0)) {
    stop("number of workers should be an integer")
  } else {
    if (strategy %in% c("multisession", "multicore", "multiprocess")) {
      message(paste0(
        "Running ",
        strategy,
        " data loading with ",
        workers,
        " workers"
      ))
      strategy <- eval(parse(text = paste0("future::", strategy)))
      future::plan(strategy, workers = workers)
    } else {
      stop("incorrect parallel computing strategy")
    }
  }
}
