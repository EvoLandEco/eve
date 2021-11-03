#' @name pdd_simulation_plot
#' @title Plotting results of a replicated pdd simulation
#' @description Function to automatically produce several plots from the results
#' of a replicated pdd simulation
#' @param result a list of extracted results
#' @param pars the parameters used in simulation being plotted
#' @param age the total simulation time
#' @param model the model used in the simulation
#' @param offset the offset method for Phi used in a pdd simulation
#' @return a plot grid generated by cow_plot
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export pdd_simulation_plot
pdd_simulation_plot <- function(result, pars, age, model, offset) {
  plotN <-
    ggplot2::ggplot(result[[2]], aes(time, N, group = rep, color = rep)) + ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") + xlab("age")

  if (model == "dsce1" | model == "dsde1") {
    plotPhi <-
      ggplot2::ggplot(result[[2]], aes(time, Phi, group = rep, color = rep)) + ggplot2::geom_line() +
      ggplot2::theme(legend.position = "none") + xlab("age") + ggplot2::geom_hline(yintercept = pars[3])
    title <-
      cowplot::ggdraw() + cowplot::draw_label(
        paste(
          model,
          ", ",
          "lambda = ",
          pars[1],
          ", ",
          "mu = ",
          pars[2],
          ", ",
          "K = ",
          pars[3],
          ", ",
          "age = ",
          age,
          ", ",
          "offset = ",
          offset
        ),
        fontface = 'bold'
      )
  } else if (model == "dsde2") {
    plotPhi <-
      ggplot2::ggplot(result[[2]], aes(time, Phi, group = rep, color = rep)) + ggplot2::geom_line() +
      ggplot2::theme(legend.position = "none") + xlab("age")
    title <-
      cowplot::ggdraw() + cowplot::draw_label(
        paste(
          model,
          ", ",
          "lambda = ",
          pars[1],
          ", ",
          "mu = ",
          pars[2],
          ", ",
          "betaN = ",
          pars[3],
          ", ",
          "betaPhi = ",
          pars[4],
          ", ",
          "gammaN = ",
          pars[5],
          ", ",
          "gammaPhi = ",
          pars[6],
          ", ",
          "age = ",
          age,
          ", ",
          "offset = ",
          offset
        ),
        fontface = 'bold'
      )
  } else if (model == "dsce2") {
    plotPhi <-
      ggplot2::ggplot(result[[2]], aes(time, Phi, group = rep, color = rep)) + ggplot2::geom_line() +
      ggplot2::theme(legend.position = "none") + xlab("age")
    title <-
      cowplot::ggdraw() + cowplot::draw_label(
        paste(
          model,
          ", ",
          "lambda = ",
          pars[1],
          ", ",
          "mu = ",
          pars[2],
          ", ",
          "betaN = ",
          pars[3],
          ", ",
          "betaPhi = ",
          pars[4],
          ", ",
          "age = ",
          age,
          ", ",
          "offset = ",
          offset
        ),
        fontface = 'bold'
      )
  }

  plotla <-
    ggplot2::ggplot(result[[2]], aes(time, lambda, group = rep, color = rep)) + ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") + xlab("age")
  plotmu <-
    ggplot2::ggplot(result[[2]], aes(time, mu, group = rep, color = rep)) + ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") + xlab("age")
  plotall <- cowplot::plot_grid(plotN, plotPhi, plotla, plotmu)

  plotwithtitle <-
    cowplot::plot_grid(title,
                       plotall,
                       ncol = 1,
                       rel_heights = c(0.1, 1))

  return(plotwithtitle)
}

#' @name edd_simulation_plot
#' @title Plotting results of a replicated edd simulation
#' @description Function to automatically produce several plots from the results
#' of a replicated edd simulation
#' @param result a list of extracted results
#' @param pars the parameters used in simulation being plotted
#' @param age the total simulation time
#' @param model the model used in the simulation
#' @return a plot grid generated by cow_plot
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_simulation_plot
edd_simulation_plot <-
  function(result,
           pars,
           age,
           model,
           offset) {
    plotN <-
      ggplot2::ggplot(result[[2]], aes(time, N, group = rep, color = rep)) + ggplot2::geom_line() +
      ggplot2::theme(legend.position = "none") + xlab("age")

    tree <-
      ggtree::ggtree(result$result_raw[, 1]$tas) + ggtree::geom_tiplab(size = 4)

    relas <-
      dplyr::bind_cols(result$result_raw[, 1]$LTT$time, result$result_raw[, 1]$las)
    relas <- relas %>% tibble::column_to_rownames(var = "...1")

    remus <-
      dplyr::bind_cols(result$result_raw[, 1]$LTT$time, result$result_raw[, 1]$mus)
    remus <- remus %>% tibble::column_to_rownames(var = "...1")

    plotlas <- ggtree::gheatmap(
      tree,
      t(relas),
      offset = 0.6,
      width = 0.8,
      colnames = FALSE,
      legend_title = "PSR"
    )

    plotmus <- ggtree::gheatmap(
      tree,
      t(remus),
      offset = 0.6,
      width = 0.8,
      colnames = FALSE,
      legend_title = "PER"
    )

    plot_bottom <- cowplot::plot_grid(plotlas, plotmus)

    plotall <-
      cowplot::plot_grid(plotN, plot_bottom, nrow = 2, rel_heights = c(1, 2))

    return(plotall)
  }

#' @name edd_plot
#' @title Generating plots for a replicated edd simulation
#' @description Function to automatically generate several plots from raw data
#' of a replicated edd simulation
#' @param raw_data a list of results generated by edd simulation function
#' @param tr true of false, to decide whether to truncate the plot to the range
#' between 0 and simulation age
#' @return a plot pack containing several plots
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_plot
edd_plot <- function(raw_data = NULL, save_plot = FALSE){
  if (is.null(raw_data)) stop("No data provided")
  if (length(raw_data) != 10) stop("Bad raw data")
  correct_names <- c("all_pars", "tes", "tas", "l_tables", "brts", "nltt",
                     "eds", "las", "mus", "linlists")
  if (!identical(names(raw_data), correct_names)) stop("Bad raw data")

  # plot lineages through time
  plot_nltt <- lapply(raw_data, edd_plot_nltt, save_plot = save_plot)

  # plot speciation rates
  plot_las <- edd_plot_las(raw_data)

  # plot extinction rates
  plot_mus <- edd_plot_mus(raw_data)

  # plot evolutionary distinctiveness-es
  plot_eds <- edd_plot_eds(raw_data)

  plot_pack <- list(plot_nltt = plot_nltt,
                    plot_las = plot_las,
                    plot_mus = plot_mus,
                    plot_eds = plot_eds
                    )

  return(plot_pack)
}

#' @name edd_plot_nltt
#' @title Generating nLTT plot for a replicated edd simulation
#' @description Function to generate normalized lineages through time plot from
#' raw data of a replicated edd simulation
#' @param raw_data a list of results generated by edd simulation function
#' @param save_plot true or false, to decide whether to save the plots to files
#' @return an plot object
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_plot_nltt
edd_plot_nltt <- function(raw_data = NULL,
                          drop_extinct = TRUE,
                          save_plot = FALSE,
                          ...) {
  if (drop_extinct = TRUE){
    df <- nLTT::get_nltt_values(raw_data$tes, dt = 0.01)
  } else {
    df <- nLTT::get_nltt_values(raw_data$tas, dt = 0.01)
  }

  lambda <- round(raw_data$all_pars$pars[[1]][1], digits = 3)
  mu <- round(raw_data$all_pars$pars[[1]][2], digits = 3)
  beta_n <- round(raw_data$all_pars$pars[[1]][3], digits = 3)
  beta_phi <- round(raw_data$all_pars$pars[[1]][4], digits = 3)
  gamma_n <- round(raw_data$all_pars$pars[[1]][5], digits = 3)
  gamma_phi <- round(raw_data$all_pars$pars[[1]][6], digits = 3)
  age <- raw_data$all_pars$age
  model <- levels(raw_data$all_pars$model)
  metric <- levels(raw_data$all_pars$metric)
  offset <- levels(raw_data$all_pars$offset)

  anno <- tibble(
    label = c(
      paste0("<span style='color:#505050'>",
             "&lambda;<sub>0</sub> = ", lambda, "<br>",
             "&mu;<sub>0</sub> = ", mu, "<br>",
             "&beta;<sub>*N*</sub> = ", beta_n, "<br>",
             "&beta;<sub>*&Phi;*</sub> = ", beta_phi, "<br>",
             "&gamma;<sub>*N*</sub> = ", gamma_n, "<br>",
             "&gamma;<sub>*&Phi;*</sub> = ", gamma_n, "</span>"),
      paste0("<span style='color:#505050'>",
             "Age = ", age, "<br>",
             "Model = ", model, "<br>",
             "Metric = ", metric, "<br>",
             "Offset = ", offset, "</span>")
    ),
    x = c(0, 0),
    y = c(0.33, .74),
    hjust = c(0, 0),
    vjust = c(0, 0),
    angle = c(0, 0)
  )

  plot_nltt <- ggplot() + geom_point(data = df, aes(t, nltt, color = id), size = I(0.1)) +
    stat_summary(data = df, aes(t, nltt),
                 fun.data = "mean_cl_boot",
                 geom = "smooth") +
    ggtitle("Average nLTT plot of phylogenies") +
    labs(x = "Normalized time", y = "Normalized number of lineages") +
    #scale_colour_ggthemr_d() +
    geom_richtext(data = anno, aes(
      x,
      y,
      label = label,
      angle = angle,
      hjust = hjust,
      vjust = vjust
    ), fill = "#E8CB9C") +
    viridis::scale_colour_viridis(discrete = TRUE, option = "A") +
    xlim(0, 1) +
    ylim(0, 1) +
    theme(legend.position = "none",
          aspect.ratio = 3 / 4)

  if (save_plot == TRUE) {
    ggsave(paste0("result/plot/nltt/",
                  lambda,"_",
                  mu,"_",
                  beta_n,"_",
                  beta_phi,"_",
                  gamma_n,"_",
                  gamma_phi,"_",
                  age,"_",
                  model,"_",
                  metric,"_",
                  offset,
                  ".png"), device = "png", dpi = "retina")

  }

  return(plot_nltt)
}



#' @name edd_plot_ltt
#' @title Generating LTT plot for a replicated edd simulation
#' @description Function to generate lineages through time plot from
#' raw data of a replicated edd simulation
#' @param raw_data a list of results generated by edd simulation function
#' @param save_plot true or false, to decide whether to save the plots to files
#' @return an plot object
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_plot_ltt
edd_plot_ltt <- function(raw_data = NULL, save_plot = FALSE, ...){
  nrep <- length(raw_data$tes)
  ltt <- raw_data$nltt
  lambda <- round(raw_data$all_pars$pars[[1]][1], digits = 3)
  mu <- round(raw_data$all_pars$pars[[1]][2], digits = 3)
  beta_n <- round(raw_data$all_pars$pars[[1]][3], digits = 3)
  beta_phi <- round(raw_data$all_pars$pars[[1]][4], digits = 3)
  gamma_n <- round(raw_data$all_pars$pars[[1]][5], digits = 3)
  gamma_phi <- round(raw_data$all_pars$pars[[1]][6], digits = 3)
  age <- raw_data$all_pars$age
  model <- levels(raw_data$all_pars$model)
  metric <- levels(raw_data$all_pars$metric)
  offset <- levels(raw_data$all_pars$offset)

  anno <- tibble(
    label = c(
      paste0("<span style='color:#505050'>",
             "&lambda;<sub>0</sub> = ", lambda, "<br>",
             "&mu;<sub>0</sub> = ", mu, "<br>",
             "&beta;<sub>*N*</sub> = ", beta_n, "<br>",
             "&beta;<sub>*&Phi;*</sub> = ", beta_phi, "<br>",
             "&gamma;<sub>*N*</sub> = ", gamma_n, "<br>",
             "&gamma;<sub>*&Phi;*</sub> = ", gamma_n, "</span>"),
      paste0("<span style='color:#505050'>",
             "Age = ", age, "<br>",
             "Model = ", model, "<br>",
             "Metric = ", metric, "<br>",
             "Offset = ", offset, "</span>")
    ),
    x = c(0, 0),
    y = c(.33, .75),
    hjust = c(0, 0),
    vjust = c(0, 0),
    angle = c(0, 0)
  )

  if (nrep != 1){
    ltt <- eve::bind_raw(ltt, nrep)
    ltt <- dplyr::bind_rows(ltt)
  }

  plot_options = TRUE

  if (plot_options == TRUE) {
    plot_ltt <-
      ggplot2::ggplot(ltt,
                      ggplot2::aes(time, num, group = as.factor(nrep), color = as.factor(nrep))) +
      ggplot2::geom_line() + ggplot2::coord_cartesian(xlim = c(0, age)) +
      ggtitle("LTT plot of phylogenies") +
      ggplot2::theme(legend.position = "none",
                     aspect.ratio = 3 / 4) +
      # geom_richtext(data = anno, aes(
      #   x,
      #   y,
      #   label = label,
      #   angle = angle,
      #   hjust = hjust,
      #   vjust = vjust
      # ), fill = "#E8CB9C") +
      viridis::scale_colour_viridis(discrete = TRUE, option = "A") +
      ggplot2::xlab("Time") + ggplot2::ylab("Number of lineages")
  } else{
    plot_ltt <-
      ggplot2::ggplot(ltt,
                      ggplot2::aes(time, num, group = as.factor(nrep), color = as.factor(nrep))) +
      ggplot2::geom_line() +
      ggtitle("LTT plot of phylogenies")
      ggplot2::theme(legend.position = "none",
                     aspect.ratio = 3 / 4) +
      viridis::scale_colour_viridis(discrete = TRUE, option = "A") +
      ggplot2::xlab("Time") + ggplot2::ylab("Number of lineages")
  }

  ggsave(paste0("result/plot/ltt/",
                lambda,"_",
                mu,"_",
                beta_n,"_",
                beta_phi,"_",
                gamma_n,"_",
                gamma_phi,"_",
                age,"_",
                model,"_",
                metric,"_",
                offset,
                ".png"), device = "png", dpi = "retina")

  return(plot_ltt)
}



#' @name edd_plot_las
#' @title Generating speciation rates plot for a replicated edd simulation
#' @description Function to generate a plot showing the transition of speciation
#' rates from raw data of a replicated edd simulation
#' @param raw_data a list of results generated by edd simulation function
#' @return an plot containing  plot
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export edd_plot_las
edd_plot_las <- function(raw_data = NULL, rep_id = 1){
  tree <-
    ggtree::ggtree(raw_data$tas[[rep_id]]) + ggtree::geom_tiplab(align = TRUE) +
    ggtree::theme_tree2()

  dat <- eve::match_raw(raw_data$las[[rep_id]],raw_data$linlists[[rep_id]])
  dat <- dplyr::bind_rows(dat)
  dat <- dplyr::bind_cols(dplyr::select(raw_data$nltt[[rep_id]], time), dat)
  dat <- dat %>%tibble::column_to_rownames(var = "time")

  gheatmap(tree, t(dat))

  return(plot_las)
}
