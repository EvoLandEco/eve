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
           device = "png",
           width = 5,
           height = 4,
           dpi = "retina") {
    if (is.null(which)) {
      stop("Plot type not specified")
    }

    save_path <- paste0("plot/", which, "/")

    check_folder(save_path, verbose = FALSE)

    if (pars_list$model == "dsce2") {
      ggplot2::ggsave(
        paste0(
          "result/",
          save_path,
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
