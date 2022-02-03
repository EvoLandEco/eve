#' extract_which
#'
#' @author Tianjian Qin
#' @return
extract_which <-
  function(raw_data = NULL,
           nrep = NULL,
           which = NULL,
           nlist = 10) {
    return(raw_data[seq(which, nlist * nrep, by = nlist)])
  }



#' bind_replication
#'
#' @author Tianjian Qin
#' @return
bind_replication <- function(raw_data = NULL, nrep = NULL) {
  if (nrep == 1)
    stop("Simulation is not replicated")

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
match_raw <- function(x = NULL, y = NULL) {
  purrr::modify2(x, y, ~ if (is_list(.x))
    fun(.x, .y)
    else
      purrr::set_names(.x, paste0('t', abs(.y))))
}



#' bind_raw
#'
#' @author Tianjian Qin
#' @return
bind_raw <- function(raw_data) {
  purrr::map(raw_data,  ~ lapply(., dplyr::bind_rows))
}



#' match_which
#'
#' @author Tianjian Qin
#' @return
match_which <- function(raw_data = NULL, which = NULL) {
  stopifnot(is.character(which))

  progress_match <-
    progressr::progressor(steps = length(raw_data$las))

  purrr::map2(
    .x = eval(parse(text = paste0(
      "raw_data$", which
    ))),
    .y = raw_data$linlists,
    .f = function(x, y) {
      progress_match()
      match_raw(x, y)
    }
  )
}


#' edd_load
#'
#' @author Tianjian Qin
#' @return
edd_load <-
  function(raw_data = NULL,
           strategy = "sequential",
           workers = NULL) {
    progressr::handlers(list(
      progressr::handler_progress(
        format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
        width    = 60,
        complete = "+"
      )
    ))

    if (is.null(workers) |
        identical(strategy, "sequential")) {
      message("Running sequential data extraction")
      message(paste0("Size of parameter sets is: ", length(raw_data)))
      message(paste0(
        "Number of replications for each parameter set is: ",
        length(raw_data$`1`$las)
      ))

      message(paste0("Matching historical states of speciation rate per lineage"))
      las <- progressr::with_progress({
        purrr::map(.x = raw_data,
                   .f = match_which,
                   which = "las")
      })

      message(paste0("Matching historical states of extinction rate per lineage"))
      mus <- progressr::with_progress({
        purrr::map(.x = raw_data,
                   .f = match_which,
                   which = "mus")
      })

      message(paste0(
        "Matching historical states of evolutionary distinctiveness per lineage"
      ))
      eds <- progressr::with_progress({
        purrr::map(.x = raw_data,
                   .f = match_which,
                   which = "eds")
      })
    } else if (!(workers %% 1 == 0)) {
      stop("number of workers should be an integer")
    } else {
      if (strategy %in% c("multisession", "multicore", "multiprocess")) {
        message(paste0("Running ",
                       strategy,
                       " loading with ",
                       workers,
                       " workers"))
        message(paste0("Size of parameter sets is: ", length(raw_data)))
        message(paste0(
          "Number of replications for each parameter set is: ",
          length(raw_data$`1`$las)
        ))

        strategy <- eval(parse(text = paste0("future::", strategy)))
        future::plan(strategy, workers = workers)

        message(paste0("Matching historical states of speciation rate per lineage"))
        las <- progressr::with_progress({
          furrr::future_map(.x = raw_data,
                            .f = match_which,
                            which = "las")
        })

        message(paste0("Matching historical states of extinction rate per lineage"))
        mus <- progressr::with_progress({
          furrr::future_map(.x = raw_data,
                            .f = match_which,
                            which = "mus")
        })

        message(
          paste0(
            "Matching historical states of evolutionary distinctiveness per lineage"
          )
        )
        eds <- progressr::with_progress({
          furrr::future_map(.x = raw_data,
                            .f = match_which,
                            which = "eds")
        })
      } else {
        stop("incorrect parallel computing strategy")
      }
    }

    # Historical states
    message("Binding data")
    hs <- lapply(list(las = las, mus = mus, eds = eds), bind_raw)

    return(hs)
  }



#' edd_merge
#'
#' @author Tianjian Qin
#' @return
edd_merge <- function(name = NULL) {
  folder_path <- paste0("result/", name)
  files <- list.files(folder_path)
  files_ordered <- gtools::mixedsort(files)
  data_path <- paste0(folder_path, "/", files_ordered)

  out <- lapply(data_path, function(x) {
    load(file = x)
    get("out")
  })

  out <- lapply(data_path, function(x) {
    print(x)
  })
}


#' edd_sim_rep
#'
#' @author Tianjian Qin
#' @return
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
edd_sim_batch <- function(nrep = 1000,
                          combo = NULL,
                          strategy = "sequential",
                          workers = NULL) {
  if (is.null(combo)) {
    stop("combo is not provided")
  }

  progress_sim <- progressr::progressor(steps = length(combo))

  if (is.null(workers) | strategy == "sequential") {
    message("Running sequential simulation")
    message(paste0("Size of parameter space is: ", length(combo)))
    message(paste0("Number of replications for each parameter set is: ", nrep))
    purrr::map(
      .x = combo,
      .f = function(x, ...) {
        progress_sim()
        edd_sim_rep(combo = x, nrep = nrep)
      },
      nrep = nrep
    )
  } else if (!(workers %% 1 == 0)) {
    stop("number of workers should be an integer")
  } else {
    if (strategy %in% c("multisession", "multicore", "multiprocess")) {
      message(paste0(
        "Running ",
        strategy,
        " simulation with ",
        workers,
        " workers"
      ))
      message(paste0("Size of parameter space is: ", length(combo)))
      message(paste0("Number of replications for each parameter set is: ", nrep))

      strategy <- eval(parse(text = paste0("future::", strategy)))
      future::plan(strategy, workers = workers)
      future_opts <- furrr::furrr_options(seed = TRUE)
      furrr::future_map(
        .x = combo,
        .f = function(x, ...) {
          progress_sim()
          edd_sim_rep(combo = x, nrep = nrep)
        },
        .options = future_opts,
        nrep = nrep
      )
    } else {
      stop("incorrect parallel computing strategy")
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
#' @param seed set random seed
#' @param strategy determine if the simulation is sequential or multi-sessioned
#' @param workers determine how many sessions are participated in the simulation
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_go
edd_go <- function(combo = NULL,
                   nrep = 1000,
                   name = NULL,
                   seed = NULL,
                   strategy = "sequential",
                   workers = NULL) {
  if (!is.null(name)) {
    if (name != "no_save") {
      check_folder(name)
    }
  } else {
    folder_name <- paste0("sim_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    check_folder(folder_name)
  }

  if (!is.null(seed) & (seed %% 1 == 0)) {
    set.seed(seed)
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
    combo$pars <- purrr::pmap(unname(combo[, 1:6]), c)
    combo <- combo[, -(1:6)]
    combo <- split(combo, seq(nrow(combo)))
    return(combo)
  } else {
    readr::write_csv2(combo, paste0("result/combo_", format(Sys.time(), "%Y%m%d_%H%M%S.csv")))
  }
}
