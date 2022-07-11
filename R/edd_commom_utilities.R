extract_parameters <- function(raw_data, unlist = FALSE) {
  if (raw_data$all_pars$model == "dsce2") {
    lambda <- round(raw_data$all_pars$pars[[1]][1], digits = 3)
    mu <- round(raw_data$all_pars$pars[[1]][2], digits = 3)
    beta_n <- round(raw_data$all_pars$pars[[1]][3], digits = 5)
    beta_phi <- round(raw_data$all_pars$pars[[1]][4], digits = 5)
    age <- raw_data$all_pars$age
    model <- raw_data$all_pars$model
    metric <- raw_data$all_pars$metric
    offset <- raw_data$all_pars$offset

    pars_list <- list(
      lambda = lambda,
      mu = mu,
      beta_n = beta_n,
      beta_phi = beta_phi,
      age = age,
      model = model,
      metric = metric,
      offset = offset
    )

    if(unlist == TRUE) {
      pars_list <- unlist(pars_list)
    }

    return(pars_list)
  } else if (raw_data$all_pars$model == "dsde2") {
    lambda <- round(raw_data$all_pars$pars[[1]][1], digits = 3)
    mu <- round(raw_data$all_pars$pars[[1]][2], digits = 3)
    beta_n <- round(raw_data$all_pars$pars[[1]][3], digits = 5)
    beta_phi <- round(raw_data$all_pars$pars[[1]][4], digits = 5)
    gamma_n <- round(raw_data$all_pars$pars[[1]][5], digits = 5)
    gamma_phi <- round(raw_data$all_pars$pars[[1]][6], digits = 5)
    age <- raw_data$all_pars$age
    model <- raw_data$all_pars$model
    metric <- raw_data$all_pars$metric
    offset <- raw_data$all_pars$offset

    pars_list <- list(
      lambda = lambda,
      mu = mu,
      beta_n = beta_n,
      beta_phi = beta_phi,
      gamma_n = gamma_n,
      gamma_phi = gamma_phi,
      age = age,
      model = model,
      metric = metric,
      offset = offset
    )

    if(unlist == TRUE) {
      pars_list <- unlist(pars_list)
    }

    return(pars_list)
  } else {
    stop("No such model")
  }
}



check_raw_data <- function(raw_data) {
  if (is.null(raw_data)) {
    stop("No raw data provided")
  }

  if (is.null(raw_data)) stop("No data provided")

  lengths <- lapply(raw_data, function(x) length(x) != 8)

  if (TRUE %in% lengths) stop("Bad raw data")

  correct_names <- c("las", "mus", "eds", "all_pars", "tes", "tas", "l_tables",
                     "ltt")

  list_names <- lapply(raw_data, function(x) !identical(names(x), correct_names))

  if (TRUE %in% list_names) {
    stop("Invalid raw data, did you forget to set history = TRUE?")
  }
}