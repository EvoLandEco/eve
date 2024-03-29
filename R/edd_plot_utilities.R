create_annotation <- function(pars_list = NULL,
                                x = c(-Inf, -Inf),
                                y = c(Inf, Inf),
                                hjust = c(0, 0),
                                vjust = c(1, 2.8),
                                angle = c(0, 0)) {
  if (pars_list$model == "dsce2") {
    anno <- tibble::tibble(
      label = c(
        paste0(
          "<span style='color:#505050'>",
          "&lambda;<sub>0</sub> = ",
          pars_list$lambda,
          "<br>",
          "&mu;<sub>0</sub> = ",
          pars_list$mu,
          "<br>",
          "&beta;<sub>*N*</sub> = ",
          pars_list$beta_n,
          "<br>",
          "&beta;<sub>*&Phi;*</sub> = ",
          pars_list$beta_phi
        ),
        paste0(
          "<span style='color:#505050'>",
          "Age = ",
          pars_list$age,
          "<br>",
          "Model = ",
          pars_list$model,
          "<br>",
          "Metric = ",
          pars_list$metric,
          "<br>",
          "Offset = ",
          pars_list$offset,
          "</span>"
        )
      ),
      x = x,
      y = y,
      hjust = hjust,
      vjust = vjust,
      angle = angle
    )

    return(anno)
  } else if (pars_list$model == "dsde2") {
    anno <- tibble::tibble(
      label = c(
        paste0(
          "<span style='color:#505050'>",
          "&lambda;<sub>0</sub> = ",
          pars_list$lambda,
          "<br>",
          "&mu;<sub>0</sub> = ",
          pars_list$mu,
          "<br>",
          "&beta;<sub>*N*</sub> = ",
          pars_list$beta_n,
          "<br>",
          "&beta;<sub>*&Phi;*</sub> = ",
          pars_list$beta_phi,
          "<br>",
          "&gamma;<sub>*N*</sub> = ",
          pars_list$gamma_n,
          "<br>",
          "&gamma;<sub>*&Phi;*</sub> = ",
          pars_list$gamma_n,
          "</span>"
        ),
        paste0(
          "<span style='color:#505050'>",
          "Age = ",
          pars_list$age,
          "<br>",
          "Model = ",
          pars_list$model,
          "<br>",
          "Metric = ",
          pars_list$metric,
          "<br>",
          "Offset = ",
          pars_list$offset,
          "</span>"
        )
      ),
      x = x,
      y = y,
      hjust = hjust,
      vjust = vjust,
      angle = angle
    )

    return(anno)
  } else {
    stop("No such model")
  }

}



save_with_parameters <-
  function(pars_list = NULL,
           plot = NULL,
           which = NULL,
           path = stop("Path not specified"),
           device = "png",
           width = 5,
           height = 4,
           dpi = "retina") {
    if (is.null(which)) {
      stop("Plot type not specified")
    }

    save_path <- file.path(path, "plot", which)

    check_path(save_path, verbose = FALSE)

    if (pars_list$model == "dsce2") {
      ggplot2::ggsave(
        filename = paste0(
          save_path,
          "/",
          which,
          "_",
          pars_list$lambda,
          "_",
          pars_list$mu,
          "_",
          pars_list$beta_n,
          "_",
          pars_list$beta_phi,
          "_",
          pars_list$age,
          "_",
          pars_list$model,
          "_",
          pars_list$metric,
          "_",
          pars_list$offset,
          ".png"
        ),
        plot = plot,
        device = device,
        width = width,
        height = height,
        dpi = dpi
      )
    } else if (pars_list$model == "dsde2") {
      ggplot2::ggsave(
        paste0(
          save_path,
          "/",
          which,
          "_",
          pars_list$lambda,
          "_",
          pars_list$mu,
          "_",
          pars_list$beta_n,
          "_",
          pars_list$beta_phi,
          "_",
          pars_list$gamma_n,
          "_",
          pars_list$gamma_phi,
          "_",
          pars_list$age,
          "_",
          pars_list$model,
          "_",
          pars_list$metric,
          "_",
          pars_list$offset,
          ".png"
        ),
        plot = plot,
        device = device,
        width = width,
        height = height,
        dpi = dpi
      )
    } else {
      stop("No such model")
    }
  }



save_with_rates_and_index_name <-
  function(rates = NULL,
           plot = NULL,
           which = NULL,
           path = stop("Path not specified"),
           device = "png",
           width = 5,
           height = 4,
           dpi = "retina") {
    if (is.null(which)) {
      stop("Plot type not specified")
    }

    save_path <- file.path(path, "plot", which)

    check_path(save_path, verbose = FALSE)

    ggplot2::ggsave(
      filename = paste0(
        save_path,
        "/",
        rates[1],
        "_",
        rates[2],
        "_",
        rates[3],
        ".png"
      ),
      plot = plot,
      device = device,
      width = width,
      height = height,
      dpi = dpi
    )
  }



save_with_rates <-
  function(rates = NULL,
           plot = NULL,
           which = NULL,
           path = stop("Path not specified"),
           device = "png",
           width = 5,
           height = 4,
           dpi = "retina") {
    if (is.null(which)) {
      stop("Plot type not specified")
    }

    save_path <- file.path(path, "plot", which)

    check_path(save_path, verbose = FALSE)

    ggplot2::ggsave(
      paste0(
        save_path,
        "/",
        which,
        "_",
        rates[1],
        ".png"
      ),
      plot = plot,
      device = device,
      width = width,
      height = height,
      dpi = dpi
    )
  }



save_with_rates2 <-
  function(rates = NULL,
           plot = NULL,
           which = NULL,
           path = stop("Path not specified"),
           device = "png",
           width = 5,
           height = 4,
           dpi = "retina") {
    if (is.null(which)) {
      stop("Plot type not specified")
    }

    save_path <- file.path(path, "plot", which)

    check_path(save_path, verbose = FALSE)

    ggplot2::ggsave(
      paste0(
        save_path,
        "/",
        which,
        "_",
        rates[1],
        "_",
        rates[2],
        "_",
        rates[3],
        "_",
        rates[4],
        ".png"
      ),
      plot = plot,
      device = device,
      width = width,
      height = height,
      dpi = dpi
    )
  }




save_with_rates_offset <-
  function(rates = NULL,
           offset = NULL,
           plot = NULL,
           which = NULL,
           path = stop("Path not specified"),
           device = "png",
           width = 5,
           height = 4,
           dpi = "retina") {
    if (is.null(which)) {
      stop("Plot type not specified")
    }

    save_path <- file.path(path, "plot", which)

    check_path(save_path, verbose = FALSE)

    ggplot2::ggsave(
      paste0(
        save_path,
        "/",
        which,
        "_",
        rates[1],
        "_",
        rates[2],
        "_",
        offset,
        ".png"
      ),
      plot = plot,
      device = device,
      width = width,
      height = height,
      dpi = dpi
    )
  }



tally_by_group <- function(data = NULL, group = NULL) {
  if (group=="metric") {
    counts <- data %>% group_by(metric, offset) %>% tally()
    return(list(rows = counts$n[1], groups = nrow(counts)))
  }
}



create_indexes_by_group <- function(tally = NULL, unlist = FALSE) {
  if (unlist) {
    indexes <- with(tally, rep(sequence(rows), each = groups) +
      (rep(rows, each = groups) * (seq(groups) -1)))
    return(indexes)
  } else {
    bases <- seq_len(tally$rows)
    indexes_list <- lapply(bases, FUN = function(x) {seq(from = x, to = tally$groups * tally$rows, by = tally$rows)})
    return(indexes_list)
  }
}



extract_tree <- function(raw_data = NULL, pars_id = NULL, rep_id = 1, drop_extinct = TRUE) {
  if (rep_id <= 0) {
    stop("Replicate ID must be positive")
  }

  if (pars_id <= 0) {
    stop("Parameter ID must be positive")
  }

  if (rep_id > length(raw_data$data$`1`$las)) {
    stop("Replicate ID exceeded number of replicates")
  }

  if (pars_id > nrow(raw_data$params)) {
    stop("Parameter ID exceeded number of parameter sets")
  }

  if (drop_extinct == TRUE) {
    tree <- raw_data$data[[pars_id]]$tes[[rep_id]]
  } else {
    tree <- raw_data$data[[pars_id]]$tas[[rep_id]]
  }

  return(tree)
}



# params is the parameter sets contained in the raw data, ... should include all the parameters that are to be extracted
find_pars_id <- function(params, pars_list = NULL, ...) {
  if (is.null(pars_list)) {
    pars_list <- list(...)
  }

  check_pars_list(pars_list)

  if (pars_list$model == "dsce2") {
    pars_id <- which(
      params$lambda == pars_list$lambda &
        params$mu == pars_list$mu &
        params$beta_n == pars_list$beta_n &
        params$beta_phi == pars_list$beta_phi &
        params$age == pars_list$age &
        params$model == pars_list$model &
        params$metric == pars_list$metric &
        params$offset == pars_list$offset
    )
  } else if (pars_list$model == "dsde2") {
    pars_id <- which(
      params$lambda == pars_list$lambda &
        params$mu == pars_list$mu &
        params$beta_n == pars_list$beta_n &
        params$beta_phi == pars_list$beta_phi &
        params$gamma_n == pars_list$gamma_n &
        params$gamma_phi == pars_list$gamma_phi &
        params$age == pars_list$age &
        params$model == pars_list$model &
        params$metric == pars_list$metric &
        params$offset == pars_list$offset
    )
  } else {
    stop("No such model")
  }

  if (length(pars_id) == 0) {
    stop("No such parameter set")
  } else if (length(pars_id) > 1) {
    stop("More than one parameter set")
  }

  return(pars_id)
}



check_pars_list <- function(pars_list) {
  if (is.null(pars_list$model)) {
    stop("Model not specified")
  }

  if (pars_list$model == "dsce2") {
    if (is.null(pars_list$lambda) |
      is.null(pars_list$mu) |
      is.null(pars_list$beta_n) |
      is.null(pars_list$beta_phi) |
      is.null(pars_list$age) |
      is.null(pars_list$metric) |
      is.null(pars_list$offset)) {
      stop("Parameter set incomplete")
    }
  } else if (pars_list$model == "dsde2") {
    if (is.null(pars_list$lambda) |
      is.null(pars_list$mu) |
      is.null(pars_list$beta_n) |
      is.null(pars_list$beta_phi) |
      is.null(pars_list$gamma_n) |
      is.null(pars_list$gamma_phi) |
      is.null(pars_list$age) |
      is.null(pars_list$metric) |
      is.null(pars_list$offset)) {
      stop("Parameter set incomplete")
    }
  } else {
    stop("No such model")
  }
}


check_pars_list <- function(pars_list) {
  if (is.null(pars_list$model)) {
    stop("Model not specified")
  }

  if (pars_list$model == "dsce2") {
    if (is.null(pars_list$lambda) |
      is.null(pars_list$mu) |
      is.null(pars_list$beta_n) |
      is.null(pars_list$beta_phi) |
      is.null(pars_list$age) |
      is.null(pars_list$metric) |
      is.null(pars_list$offset)) {
      stop("Parameter set incomplete")
    }
  } else if (pars_list$model == "dsde2") {
    if (is.null(pars_list$lambda) |
      is.null(pars_list$mu) |
      is.null(pars_list$beta_n) |
      is.null(pars_list$beta_phi) |
      is.null(pars_list$gamma_n) |
      is.null(pars_list$gamma_phi) |
      is.null(pars_list$age) |
      is.null(pars_list$metric) |
      is.null(pars_list$offset)) {
      stop("Parameter set incomplete")
    }
  } else {
    stop("No such model")
  }
}


# quickly find the replicate ID of a row in the statistics table
find_rep_id <- function(row_id, nrep) {
  if (row_id <= 0) {
      stop("Row ID must be positive")
  }

  if (nrep <= 0) {
      stop("Number of replicates must be positive")
  }

  rep_id <- row_id %% nrep

  return(rep_id)
}



pars_to_title <- function(pars_list) {
  pars <- unlist(pars_list$pars)
  if (pars_list$model == "dsce2") {
    lambda <- as.character(pars[1])
    mu <- as.character(pars[2])
    beta_n <- as.character(pars[3])
    beta_phi <- as.character(pars[4])

    plot_title <- bquote(lambda[italic(0)] ~ "=" ~ .(lambda) ~ mu[italic(0)] ~ "=" ~ .(mu) ~ beta[italic(N)] ~ "=" ~ .(beta_n) ~ beta[italic(Phi)] ~ "=" ~ .(beta_phi))
  } else if (pars_list$model == "dsde2") {
    lambda <- as.character(pars[1])
    mu <- as.character(pars[2])
    beta_n <- as.character(pars[3])
    beta_phi <- as.character(pars[4])
    gamma_n <- as.character(pars[5])
    gamma_phi <- as.character(pars[6])

    plot_title <- bquote(lambda[italic(0)] ~ "=" ~ .(lambda) ~ mu[italic(0)] ~ "=" ~ .(mu) ~ beta[italic(N)] ~ "=" ~ .(beta_n) ~ beta[italic(Phi)] ~ "=" ~ .(beta_phi) ~ gamma[italic(N)] ~ "=" ~ .(gamma_n) ~ gamma[italic(Phi)] ~ "=" ~ .(gamma_phi))
  } else {
    stop("No such model")
  }

  return(plot_title)
}


# Use with params argument in edd_stat functions
pars_to_title2 <- function(params) {
  if (params$model == "dsce2") {
    plot_title <- bquote(lambda[italic(0)] ~ "=" ~ .(params$lambda) ~ mu[italic(0)] ~ "=" ~ .(params$mu) ~ beta[italic(N)] ~ "=" ~ .(params$beta_n) ~ beta[italic(Phi)] ~ "=" ~ .(params$beta_phi))
  } else if (params$model == "dsde2") {
    plot_title <- bquote(lambda[italic(0)] ~ "=" ~ .(params$lambda) ~ mu[italic(0)] ~ "=" ~ .(params$mu) ~ beta[italic(N)] ~ "=" ~ .(params$beta_n) ~ beta[italic(Phi)] ~ "=" ~ .(params$beta_phi) ~ gamma[italic(N)] ~ "=" ~ .(params$gamma_n) ~ gamma[italic(Phi)] ~ "=" ~ .(params$gamma_phi))
  } else {
    stop("No such model")
  }

  return(plot_title)
}


index_name_to_title <- function(name) {
  if (name == "J_One") {
    return("J-One balance index")
  } else if (name == "Gamma") {
    return("Gamma statistic")
  } else if (name == "PD") {
    return("Phylogenetic diversity")
  } else if (name == "MBL") {
    return("Mean branch length")
  } else if (name == "MNTD") {
    return("Mean nearest taxon distance")
  } else if (name == "ERE") {
    return("Speciation rate evenness")
  } else if (name == "MPD") {
    return("Mean pairwise distance")
  } else if (name == "SR") {
    return("Number of lineages")
  } else if (name == "Rogers") {
    return("Rogers J index of (im)balance")
  } else if (name == "Colless") {
    return("Colless index of (im)balance")
  } else if (name == "Sackin") {
    return("Sackin index of (im)balance")
  } else if (name == "TCI") {
    return("Total cophenetic index")
  } else if (name == "B1") {
    return("B1 index of (im)balance")
  } else {
    stop("No such statistic or statistic not yet implemented")
  }
}