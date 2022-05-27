#' extract_which
#'
#' @author Tianjian Qin
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
match_raw <- function(x = NULL, y = NULL) {
  purrr::modify2(x, y, ~ if (is_list(.x))
    fun(.x, .y)
    else
      purrr::set_names(.x, paste0('t', abs(.y))))
}



#' bind_raw
#'
#' @author Tianjian Qin
bind_raw <- function(raw_data) {
  purrr::map(raw_data,  ~ lapply(., dplyr::bind_rows))
}



#' match_which
#'
#' @author Tianjian Qin
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


#' edd_load_split
#'
#' @author Tianjian Qin
edd_load_split <-
  function(raw_data = NULL,
           strategy = "sequential",
           workers = 1) {
    progressr::handlers(list(
      progressr::handler_progress(
        format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
        width    = 60,
        complete = "+"
      )
    ))

    check_parallel_arguments(workers, strategy)

    message(paste0("Size of parameter sets is: ", length(raw_data)))
    message(paste0(
      "Number of replications for each parameter set is: ",
      length(raw_data$`1`$las)
    ))

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

    message(paste0(
      "Matching historical states of evolutionary distinctiveness per lineage"
    ))
    eds <- progressr::with_progress({
      furrr::future_map(.x = raw_data,
                        .f = match_which,
                        which = "eds")
    })

    # Historical states
    message("Merging historical states")
    hs <- lapply(list(las = las, mus = mus, eds = eds), bind_raw)

    return(hs)
  }



#' edd_merge
#'
#' @author Tianjian Qin
edd_merge <- function(name = NULL) {
  folder_path <- paste0("result/", name)
  files <- list.files(folder_path)
  files_ordered <- gtools::mixedsort(files)
  data_path <- paste0(folder_path, "/", files_ordered)

  message("Merging splitted datasets")
  out <- lapply(data_path, function(x) {
    load(file = x)
    get("out")
  })

  names(out) <- 1:length(files)

  return(out)
}



#' edd_load
#'
#' @author Tianjian Qin
#' @export edd_load
edd_load <- function(name = NULL) {
  merged_data <- edd_merge(name)
  loaded_data <- edd_load_split(merged_data)

  return(loaded_data)
}
