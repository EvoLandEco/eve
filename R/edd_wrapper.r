#' edd_sim_rep
#' @export edd_sim_rep
#' @author Tianjian Qin
edd_sim_rep <-
  function(combo = NULL, history = FALSE, verbose = FALSE, nrep = 5) {
    if (nrep < 2)
      stop("Must have more than 1 replicates")

    combo <- as.data.frame(combo, col.names = NULL)

    # do replicated simulations
    raw_data <- replicate(nrep, {
      edd_sim(pars = unlist(combo$pars),
                   age = combo$age,
                   model = combo$model,
                   metric = combo$metric,
                   offset = combo$offset,
                   history = history,
                   verbose = verbose)
    })

    if (history == TRUE) {
      out <- list(
        all_pars = combo,
        tes = raw_data[1,],
        tas = raw_data[2,],
        l_tables = raw_data[3,],
        ltt = raw_data[4,],
        eds = raw_data[5,],
        las  = raw_data[6,],
        mus  = raw_data[7,],
        linlists = raw_data[8,]
      )
    } else {
      out <- list(
        all_pars = combo,
        tes = raw_data[1,],
        tas = raw_data[2,],
        l_tables = raw_data[3,],
        ltt = raw_data[4,]
      )
    }


    return(out)
  }



#' edd_sim_batch
#'
#' @author Tianjian Qin
edd_sim_batch <- function(combo = NULL,
                          history = FALSE,
                          nrep = 1000,
                          strategy = "sequential",
                          workers = NULL,
                          verbose = FALSE) {
  if (is.null(combo)) {
    stop("combo is not provided")
  }

  progress_sim <- progressr::progressor(steps = length(combo))

  check_parallel_arguments(strategy, workers, verbose)

  message(paste0("Size of parameter space is: ", length(combo)))
  message(paste0("Number of replications for each parameter set is: ", nrep))

  future_opts <- furrr::furrr_options(seed = TRUE)
  furrr::future_map(
    .x = combo,
    .f = function(x, ...) {
      progress_sim()
      edd_sim_rep(combo = x, nrep = nrep)
    },
    .options = future_opts,
    nrep = nrep,
    history = history,
    verbose = verbose
  )
}



#' @name edd_go
#' @title Starting a edd simulation
#' @description Function to start a replicated edd simulation along with a given
#' parameter space
#' @param nrep number of replication
#' @param combo parameter space
#' @param name name of the folder to save results, if left blank, a folder will
#' be created according to current time and date. Set it to "no_save" to disable
#' result saving function, an object containing all the results will be returned
#' instead.
#' @param seed set random seed
#' @param strategy determine if the simulation is sequential or multi-sessioned
#' @param workers determine how many sessions are participated in the simulation
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_go
edd_go <- function(combo = NULL,
                   history = FALSE,
                   nrep = 1000,
                   name = NULL,
                   seed = NULL,
                   strategy = "sequential",
                   workers = NULL,
                   verbose = FALSE) {
  if (!is.null(name)) {
    if (name != "no_save") {
      check_folder(name)
    }
  } else {
    folder_name <- paste0("sim_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    check_folder(folder_name)
  }

  if (!is.null(seed)) {
    if (seed %% 1 == 0) {
      set.seed(seed)
    } else {
      stop("must provide a valid random seed")
    }
  }

  progressr::handlers(list(
    progressr::handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 60,
      complete = "+"
    )
  ))

  out <- progressr::with_progress({
    edd_sim_batch(
      combo = combo,
      history = history,
      verbose = verbose,
      nrep = nrep,
      strategy = strategy,
      workers = workers
    )
  })

  if (!is.null(name)) {
    if (name != "no_save") {
      save_result(out, name)
    } else {
      return(out)
    }
  } else {
    save_result(out, folder_name)
  }
}



#' @name edd_combo_maker
#' @title Making a parameter space for edd simulation
#' @description Function to conveniently make a parameter space for edd
#' simulation
#' @param ... parameters that form a parameter space
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_combo_maker
edd_combo_maker <- function(save_file = FALSE, ...) {
  pars <- list(...)

  if (length(names(pars)) != 8 & length(names(pars)) != 10) {
    stop("incorrect number of parameters")
  }

  if (length(pars$model) > 1) {
    stop("only one model allowed for one batch simulation run")
  }

  if (pars$model == "dsce2") {
    if (length(names(pars)) == 8) {
      pars_list <- c("la",
                     "mu",
                     "beta_n",
                     "beta_phi",
                     "age",
                     "model",
                     "metric",
                     "offset")
      pars_identical <- identical(names(pars), pars_list)
      stopifnot("incorrect parameters or incorrect order" = pars_identical)
    } else {
      stop("incorrect parameters for dsce2 model")
    }
  } else if (pars$model == "dsde2") {
    if (length(names(pars)) == 10) {
      pars_list <- c(
        "la",
        "mu",
        "beta_n",
        "beta_phi",
        "gamma_n",
        "gamma_phi",
        "age",
        "model",
        "metric",
        "offset"
      )
      pars_identical <- identical(names(pars), pars_list)
      stopifnot("incorrect parameters or incorrect order" = pars_identical)
    } else {
      stop("incorrect parameters for dsde2 model")
    }
  } else{
    stop("incorrect model")
  }

  combo <- expand.grid(...)

  if (save_file == FALSE){
    if (pars$model == "dsde2") {
      combo$pars <- purrr::pmap(unname(combo[, 1:6]), c)
      combo <- combo[, -(1:6)]
      combo <- split(combo, seq(nrow(combo)))
    } else if (pars$model == "dsce2") {
      combo$pars <- purrr::pmap(unname(combo[, 1:4]), c)
      combo <- combo[, -(1:4)]
      combo <- split(combo, seq(nrow(combo)))
    }
    return(combo)
  } else {
    readr::write_csv2(combo, paste0("result/combo_", format(Sys.time(), "%Y%m%d_%H%M%S.csv")))
  }
}
