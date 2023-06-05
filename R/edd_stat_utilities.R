#' @name calculate_CI_alt
#' @param brts_list A list containing branching times of a simulation
#' @param tt Time table
#' @param alpha Resolution
#' @title Calculating confidence interval for LTT plot (reversed time scale)
#' @description Function to calculate confidence interval for LTT plot, its an Rcpp alternative version of the
#' native C++ function, used for debugging
#' @author Thijs Janzen; Tianjian Qin
calculate_CI_alt <- Rcpp::cppFunction('
NumericMatrix get_stats_cpp(const List& brts_list,
                            const NumericVector& tt,
                            const float alpha) {
  NumericMatrix output(tt.size(), 5);

  size_t min_index = brts_list.size() * alpha/2;
  size_t max_index = brts_list.size() * (1 - alpha/2);

  for (int i = (tt.size() - 1); i >= 0; --i) {
    float focal_t = tt[i];
    std::vector<size_t> found_lin(brts_list.size());
    float mean = 0.f;
    NumericVector focal_brts;
    for (int j = 0; j < brts_list.size(); ++j) {
      focal_brts = brts_list[j]; // temp
      size_t cnt = focal_brts.size() - 1;
      while(focal_brts[cnt] > focal_t) {
        cnt--;
      }
      found_lin[j] = cnt; // two at first timepoint, crown age.
      mean += found_lin[j];
    }
    mean *= 1.0 / brts_list.size();

    std::sort(found_lin.begin(), found_lin.end());
    float median = found_lin[ found_lin.size() / 2 + 1 ];

    if (found_lin.size() % 2 == 0) {
      median = 0.5 * (found_lin[ found_lin.size() / 2] +
                      found_lin[ found_lin.size() / 2 + 1]);
    }

    float minalpha = found_lin.at( min_index );
    float maxalpha = found_lin.at( max_index );

    NumericVector add(5);
    add(0) = focal_t;
    add(1) = median;
    add(2) = minalpha;
    add(3) = maxalpha;
    add(4) = mean;

    output(i, _) = add;
  }
  return output;
}')


calculate_tree_stats <- function(phy = NULL, min_size = 3, method = "treestats", metric = "Aldous") {
  if (!ape::is.ultrametric.phylo(phy)) {
    stop("Ultrametric tree required, considering using tree with only extant taxa")
  }

  if ((phy$Nnode + 1) >= min_size) {
    if (method == "treestats") {
      if (metric == "Aldous") {
        return(treestats::beta_statistic(phy))
      } else if (metric == "Sackin") {
        return(treestats::sackin(phy, normalization = "yule"))
      } else if (metric == "Colless") {
        return(treestats::colless(phy, normalization = "yule"))
      } else if (metric == "Blum") {
        return(treestats::blum(phy))
      } else if (metric == "J-One") {
        return(treestats::j_one(phy))
      } else if (metric == "TCI") {
        return(treestats::tot_coph(phy, normalization = "yule"))
      } else if (metric == "Steps") {
        return(treestats::imbalance_steps(phy, normalize = TRUE))
      } else if (metric == "Branch") {
        return(calc_branch_colless(phy))
      } else if (metric == "B1") {
        return(treestats::b1(phy, normalization = "none"))
      } else if (metric == "B2") {
        return(treestats::b2(phy, normalization = "yule"))
      } else if (metric == "Gamma") {
        return(treestats::gamma_statistic(phy))
      } else if (metric == "MBL") {
        return(treestats::mean_branch_length(phy))
      } else if (metric == "PD") {
        return(treestats::phylogenetic_diversity(phy))
      } else if (metric == "MNTD") {
        return(treestats::mntd(phy))
        stop("Invalid metric")
      }
    } else {
      stop("No such method")
    }
  } else {
    return(NA)
  }
}


transform_data <- function(stat) {
  if (!is.null(stat$offset)) {
    stat <- dplyr::mutate(stat,
                          offset = dplyr::case_when(offset == "none" ~ "None",
                                                    offset == "simtime" ~ "Simulation time",
                                                    offset == "spcount" ~ "Species count",
                                                    offset == "both" ~ "Both"))
    stat$offset <- factor(stat$offset, levels = c("None", "Simulation time", "Species count", "Both"))
  }
  if (!is.null(stat$lambda)) {
    stat$lambda <- as.factor(stat$lambda)
  }
  if (!is.null(stat$mu)) {
    stat$mu <- as.factor(stat$mu)
  }
  if (!is.null(stat$beta_n)) {
    stat$beta_n <- as.factor(stat$beta_n)
  }
  if (!is.null(stat$beta_phi)) {
    stat$beta_phi <- as.factor(stat$beta_phi)
  }

  if (stat$model[1] == "dsde2") {
    if (!is.null(stat$gamma_n)) {
      stat$gamma_n <- as.factor(stat$gamma_n)
    }
    if (!is.null(stat$gamma_phi)) {
      stat$gamma_phi <- as.factor(stat$gamma_phi)
    }
  }

  return(stat)
}


reverse_l_table <- function(l_table, age) {
  l_table[, 1] <- age - l_table[, 1]
  return(l_table)
}


get_stats_names <- function(raw_data, stats) {
  stats_names <- colnames(stats)
  model <- raw_data$params$model[1]
  if (model == "dsce2") {
    stats_names <- stats_names[-(1:8)]
  } else if (model == "dsde2") {
    stats_names <- stats_names[-(1:10)]
  } else {
    stop("Model not supported")
  }
  return(stats_names)
}


check_stats_names <- function(raw_data, stats, name) {
  if (!all(is.character(name))) {
    stop("Statistic name must be character")
  }

  duplicates <- duplicated(name)
  if (any(duplicates)) {
    stop("Duplicated statistic names: ", paste(name[duplicates], collapse = ", "))
  }

  stats_names <- get_stats_names(raw_data, stats)

  result <- name %in% stats_names

  if (FALSE %in% result) {
    stop("Statistic(s) does not exist. Available statistic(s): ", paste(stats_names, collapse = ", "))
  }
}


find_best_rep_ids <- function(raw_data = NULL, method = "euclidean", metrics = NULL) {
  stats <- edd_stat_cached(raw_data$data)

  if (is.null(metrics)) {
    message("No metrics specified, using all")
    metrics <- get_stats_names(raw_data, stats)
  } else {
    check_stats_names(raw_data, stats, metrics)
  }

  excluded_metrics <- setdiff(get_stats_names(raw_data, stats), metrics)

  rep_ids <- stats %>%
    dplyr::mutate(group = interaction(lambda, mu, beta_n, beta_phi, age, model, metric, offset)) %>%
    dplyr::select(-dplyr::all_of(excluded_metrics)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(rep_id = dplyr::row_number()) %>%
    na.omit() %>%
    dplyr::mutate(distance = calculate_tree_distance_cached(dplyr::cur_data()[, dplyr::all_of(metrics)], method)) %>%
    dplyr::slice_min(n = 1, order_by = distance) %>%
    dplyr::ungroup() %>%
    dplyr::select(-group) %>%
    dplyr::mutate(pars_id = dplyr::row_number())

  # cols <- stats %>%
  #   dplyr::select(-lambda, -mu, -beta_n, -beta_phi, -age, -model, -metric, -offset) %>%
  #   dplyr::mutate_all(~(scale(.) %>% as.vector))
  #
  # col_means <- cols %>% dplyr::mutate_all(~ .x - mean(.x, na.rm = TRUE)) %>% dplyr::ungroup() %>% dplyr::select(-group)
  #
  # stats$distance <- apply(col_means, 1, function(x) sqrt(sum(x^2)))

  # rep_ids <- stats %>%
  #   dplyr::mutate(rep_id = dplyr::row_number()) %>%
  #   dplyr::slice_min(n = 1, order_by = distance) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(-group) %>%
  #   dplyr::mutate(pars_id = row_number())

  return(rep_ids)
}


calculate_tree_distance <- function(stats = NULL, method = "euclidean") {
  if (method == "euclidean") {
    mean_stats <- colMeans(stats)
    return(sqrt(rowSums((stats - mean_stats)^2)))
  } else if (method == "manhattan") {
    mean_stats <- colMeans(stats)
    return(rowSums(abs(stats - mean_stats)))
  } else if (method == "mahalanobis") {
    mean_stats <- colMeans(stats)
    return(mahalanobis(stats, mean_stats, cov(stats)))
  } else {
    stop("The specified method is not supported. Choose from 'euclidean', 'manhattan', 'mahalanobis'.")
  }
}