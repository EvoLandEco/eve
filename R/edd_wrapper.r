#' extract_edd_result
#'
#' @author Tianjian Qin
#' @return
#' @export
extract_edd_result <- function(result, nrep, which, nlist = 8) {
  out <-
    return(result[seq(which, nlist * nrep, by = nlist)])
}

#' bind_raw
#'
#' @author Tianjian Qin
#' @return
#' @export
bind_raw <- function(raw_data = NULL, nrep = NULL) {
  if (nrep == 1) stop("Simulation is not replicated")
  binded_data <- lapply(raw_data, as.data.frame)

  for (i in 1:nrep) {
    binded_data[[i]]["nrep"] <- rep(i, times = nrow(binded_data[[i]]))
  }

  return(binded_data)
}

#' match_raw
#'
#' @author Tianjian Qin
#' @return
#' @export
#' @author Tianjian Qin
#' @return
#' @export
match_raw <- function(x = NULL, y = NULL){
 purrr::modify2(x, y, ~if(is_list(.x)) fun(.x, .y)
                  else purrr::set_names(.x, paste0('t', abs(.y))))
}

#' edd_sim_rep
#'
#' @author Tianjian Qin
#' @return
#' @export
edd_sim_rep <-
  function(combo = NULL, nrep = 5) {
    if (nrep < 2)
      stop("Must have more than 1 replicates")

    combo <- as.data.frame(combo, col.names = NULL)

    # do replicated simulations
    raw_data <- replicate(nrep, {
      DDD::edd_sim(unlist(combo$pars),
                   combo$age,
                   combo$model,
                   combo$metric,
                   combo$offset)
    })

    out <- list(
      all_pars = combo,
      tes = raw_data[1,],
      tas = raw_data[2,],
      l_tables = raw_data[3,],
      brts = raw_data[4,],
      nltt = raw_data[5,],
      eds = raw_data[6,],
      las  = raw_data[7,],
      mus  = raw_data[8,],
      linlists = raw_data[9,]
    )

    return(out)
  }

#' edd_sim_batch
#'
#' @author Tianjian Qin
#' @return
#' @export
edd_sim_batch <- function(nrep = 1000,
                          combo = NULL,
                          strategy = future::sequential,
                          workers = NULL) {
  if (is.null(combo)) {
    stop("combo is not provided")
  }

  progress_sim <- progressr::progressor(steps = length(combo))

  if (is.null(workers) | identical(strategy, future::sequential)) {
    message("Running sequential simulation")
    message(paste0("Number of parameter sets is: ", length(combo)))
    message(paste0("Number of replication for each parameter set is: ", nrep))
    purrr::map(
      .x = combo,
      .f = function(x, ...) {
        progress_sim()
        eve::edd_sim_rep(combo = x, nrep = nrep)
      },
      nrep = nrep
    )
  } else if (!(workers %% 1 == 0)) {
    stop("number of workers should be an integer")
  } else {
    if (identical(strategy, future::multisession)) {
      message(paste0(
        "Running multisession simulation with ",
        workers,
        " workers"
      ))
      message(paste0("Number of parameter sets is: ", length(combo)))
      message(paste0("Number of replication for each parameter set is: ", nrep))
      future::plan(strategy, workers = workers)
      future_opts <- furrr::furrr_options(seed = TRUE)
      furrr::future_map(
        .x = combo,
        .f = function(x, ...) {
          progress_sim()
          eve::edd_sim_rep(combo = x, nrep = nrep)
        },
        .options = future_opts,
        nrep = nrep
      )
    } else {
      stop("incorrect simulation strategy")
    }
  }
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
#' @param strategy determine if the simulation is sequential or multi-sessioned
#' @param workers determine how many sessions are participated in the simulation
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_go
edd_go <- function(combo = NULL,
                   nrep = 1000,
                   name = NULL,
                   strategy = future::sequential,
                   workers = NULL) {
  if (!is.null(name)) {
    if (name != "no_save") {
      eve::check_folder(name)
    }
  } else {
    eve::check_folder(paste0("sim_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }

  progressr::handlers(list(
    progressr::handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 60,
      complete = "+"
    )
  ))

  out <- progressr::with_progress({
    eve::edd_sim_batch(
      combo = combo,
      nrep = nrep,
      strategy = strategy,
      workers = workers
    )
  })

  if (!is.null(name)) {
    if (name != "no_save") {
      eve::save_result(out, name)
    } else {
      return(out)
    }
  } else {
    eve::save_result(out, paste0("sim_", format(Sys.time(), "%Y%m%d_%H%M%S")))
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
    combo$pars <- purrr::pmap(unname(combo[, 1:6]), c)
    combo <- combo[, -(1:6)]
    combo <- split(combo, seq(nrow(combo)))
    return(combo)
  } else {
    write_csv2(combo, paste0("result/combo_", format(Sys.time(), "%Y%m%d_%H%M%S.csv")))
  }
}
