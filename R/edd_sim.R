edd_pars_check <- function(pars, age, model, metric, offset) {
  # pars range check
  if (pars[1] <= 0) {
    stop("speciation rate should be positive")
  }
  if (pars[2] < 0) {
    stop("extinction rate should be none-negative")
  }
  # pars and model match check
  if (model == "dsce2" && length(pars) != 4) {
    stop("this model requires four parameters")
  }
  if (model == "dsde2" && length(pars) != 6) {
    stop("this model requires six parameters")
  }
  # metric and offset match check
  if (metric != "pd" && offset != "none") {
    stop("only pd metric has offset methods")
  }
  if (age <= 0) {
    stop("age should be positive")
  }
}

edd_message_info <- function(pars, age, model, metric, offset) {
  if (model == "dsce2") {
    message("Simulation Info:")
    message(paste0("Lambda = ", pars[1]))
    message(paste0("Mu = ", pars[2]))
    message(paste0("Beta_N = ", pars[3]))
    message(paste0("Beta_Phi= ", pars[4]))
    message(paste0("Age = ", age))
    message(paste0("Model = ", model))
    message(paste0("Metric = ", metric))
    message(paste0("Offset = ", offset))
  } else if (model == "dsde2") {
    message("Simulation Info:")
    message(paste0("Lambda = ", pars[1]))
    message(paste0("Mu = ", pars[2]))
    message(paste0("Beta_N = ", pars[3]))
    message(paste0("Beta_Phi= ", pars[4]))
    message(paste0("Gamma_N = ", pars[5]))
    message(paste0("Gamma_Phi = ", pars[6]))
    message(paste0("Age = ", age))
    message(paste0("Model = ", model))
    message(paste0("Metric = ", metric))
    message(paste0("Offset = ", offset))
  }
}

edd_update_lamu <- function(ed, ed_max, params, model) {
  num <- params[1]
  la0 <- params[2]
  mu0 <- params[3]
  beta_num <- params[4]
  beta_phi <- params[5]
  if (model == "dsce2") {
    if (length(params) != 5) {
      stop("incorrect parameter(s)")
    }
    # dependent speciation, constant extinction
    if (beta_phi < 0) {
      newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed)
    } else {
      newlas <- pmin(la0 + beta_num * num + beta_phi * ed, la0 + beta_num * num + beta_phi * ed_max)
      newlas <- pmax(0, newlas)
    }
    newmus <- rep(mu0, length(newlas))
  } else if (model == "dsde2") {
    if (length(params) != 7) {
      stop("incorrect parameter(s)")
    }
    # dependent speciation, dependent extinction
    gamma_num <- params[6]
    gamma_phi <- params[7]
    if (beta_phi < 0) {
      newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed)
    } else {
      newlas <- pmin(la0 + beta_num * num + beta_phi * ed, la0 + beta_num * num + beta_phi * ed_max)
      newlas <- pmax(0, newlas)
    }
    if (gamma_phi < 0) {
      newmus <- pmax(0, mu0 + gamma_num * num + gamma_phi * ed)
    } else {
      newmus <- pmin(mu0 + gamma_num * num + gamma_phi * ed, mu0 + gamma_num * num + gamma_phi * ed_max)
      newmus <- pmax(0, newmus)
    }
  }
  return(list(newlas = newlas, newmus = newmus))
}

edd_get_edmax <-
  function(num, l_table, t, metric, offset, converter = "cpp") {
    if (metric == "ed") {
      if (converter == "cpp") {
        ed_max <- as.vector(rep(2 * t, num))
      } else {
        ed_max <- as.vector(rep(2 * t, num))
      }
    } else if (metric == "nnd") {
      ed_max <- as.vector(rep(2 * t, num))
    } else if (metric == "pd") {
      if (offset == "none") {
        if (converter == "cpp") {
          ed_max <- rep(as.vector(L2Phi_cpp(l_table, t, metric)), num)
        } else {
          ed_max <- rep(as.vector(L2Phi(l_table, t, metric)), num)
        }
      } else if (offset == "simtime") {
        if (converter == "cpp") {
          ed_max <-
            rep(as.vector(L2Phi_cpp(l_table, t, metric) - t), num)
        } else {
          ed_max <- rep(as.vector(L2Phi(l_table, t, metric) - t), num)
        }
      } else {
        stop("no such offset method")
      }
    } else {
      stop("no such metric")
    }
    return(ed_max)
  }

edd_get_ed <-
  function(num, l_table, t, metric, offset, converter = "cpp") {
    if (metric == "ed") {
      if (converter == "cpp") {
        ed <- as.vector(L2ED_cpp(l_table, t))
      } else {
        ed <- as.vector(L2ED(l_table, t))
      }
    } else if (metric == "nnd") {
      ed <- as.vector(L2NND(l_table, t))
    } else if (metric == "pd") {
      if (offset == "none") {
        if (converter == "cpp") {
          ed <- rep(as.vector(L2Phi_cpp(l_table, t, metric)), num)
        } else {
          ed <- rep(as.vector(L2Phi(l_table, t, metric)), num)
        }
      } else if (offset == "simtime") {
        if (converter == "cpp") {
          ed <- rep(as.vector(L2Phi_cpp(l_table, t, metric) - t), num)
        } else {
          ed <- rep(as.vector(L2Phi(l_table, t, metric) - t), num)
        }
      } else {
        stop("no such offset method")
      }
    } else {
      stop("no such metric")
    }
    return(ed)
  }

edd_sum_rates <- function(las, mus) {
  return(sum(las) + sum(mus))
}

edd_sample_event <- function(las, mus, linlist) {
  events <- 1:(2 * length(linlist))
  return(sample2(events, 1, prob = c(las, mus)))
}

#' @title Function to simulate the evolutionary distinctiveness dependent
#' diversification process
#' @description Simulating the evolutionary distinctiveness dependent
#' diversification process
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to
#' lambda (speciation rate) \cr \code{pars[2]} corresponds to mu (extinction
#' rate) \cr \code{pars[3]} corresponds to beta_num (coefficient for species
#' number effect on speciation) \cr \code{pars[4]} corresponds to beta_phi
#' (coefficient for evolutionary distinctiveness effect on speciation)
#' \cr \code{pars[5]} corresponds to gamma_num (coefficient for species number
#' effect on speciation) \cr \code{pars[6]} corresponds to gamma_phi
#' (coefficient for evolutionary distinctiveness effect on extinction)
#' @param age Sets the crown age for the simulation
#' @param model Sets the model of diversity-dependence: \cr \code{model ==
#' dsce2} : linear dependence in speciation rate with parameters beta_num and
#' beta-phi\cr \code{model == dsde2} : linear dependence
#' in both speciation rate and extinction rate with parameters beta_num,
#' beta_phi, gamma_num and gamma_phi
#' @param metric "pd" , "ed" or "nnd", Specifies which evolutionary distinctiveness
#' metric should be used.
#' @param offset Specifies which method to use to offset the impact of tree age
#' and the collinearity between pd and species richness. "none" for no offset
#' method; "simtime" for deducting tree age from pd value; "spcount" for
#' dividing pd value by species richness; "both" for applying both "simtime" and
#' "spcount", by firstly deducting tree age and then dividing by species richness
#' @param history Logical, indicating whether to record the historical states
#' (of the rates and ED/PD values)
#' @param verbose Logical, for debugging purpose, indicating whether to print
#' simulation info at each step in the console, and save running time to a file
#' @param converter Either "cpp" or "r", choose which version of L2phylo to use.
#' @param size_limit The maximum number of species allowed in the simulation.
#' @return \item{ out }{ A list with the following nine elements: The first
#' element is the tree of extant species in phylo format \cr The second element
#' is the tree of all species, including extinct species, in phylo format \cr
#' The third element is a matrix of all species where \cr - the first column is
#' the time at which a species is born \cr - the second column is the label of
#' the parent of the species; positive and negative values only indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values only indicate whether the species belongs to the left or
#' right crown lineage \cr - the fourth column is the time of extinction of the
#' species. If this equals -1, then the species is still extant.\cr The fourth
#' element is the set of branching times of the tree of extant species.\cr The
#' fifth element is the lineage-through-time plot. \cr The sixth element is a
#' list of all evolutionary distinctiveness values of all lineages. \cr The
#' seventh element is a list of all the speciation rates of all lineages at all
#' the time steps. \cr The eighth element is a list of all the extinction rates
#' of all lineages at all the time steps. \cr The ninth element is a list of
#' all the lineages at all the time steps.}
#' @author Tianjian Qin, Rampal S. Etienne
#' @keywords models
#' @examples
#' edd_sim(
#'   pars = c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001), age = 6,
#'   model = "dsde2", metric = "ed", offset = "none"
#' )
#' @export edd_sim
edd_sim <- function(pars,
                    age,
                    model = "dsce2",
                    metric = "ed",
                    offset = "none",
                    history = TRUE,
                    verbose = FALSE,
                    converter = "cpp",
                    size_limit = 2000) {
  edd_pars_check(pars, age, model, metric, offset)
  if (verbose == TRUE) {
    edd_message_info(pars, age, model, metric, offset)
  }

  done <- 0
  while (done == 0) {
    i <- 1
    t <- rep(0, 1)
    l_table <- matrix(0, 2, 4)
    t[1] <- 0
    num <- 2
    l_table[1, 1:4] <- c(0, 0, -1, -1)
    l_table[2, 1:4] <- c(0, -1, 2, -1)
    linlist <- c(-1, 2)
    new_lin <- 2
    params <- c(num, pars)
    ed <- c(0, 0)
    ed_max <-
      edd_get_edmax(num, l_table, age, metric, offset, converter)
    lamu <- edd_update_lamu(ed, ed_max, params, model)

    if (history == TRUE) {
      eds <- list(ed)
      las <- list(lamu$newlas)
      mus <- list(lamu$newmus)
      linlists <- list(linlist)
    }

    t[i + 1] <-
      t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))

    if (verbose == TRUE) {
      message(paste0("Simulation step ", i, " started"))
      message("All variables initialized")
      message(paste0("The first time interval is ", t[i + 1]))
      times <- rep(0, 1)
    }

    while (t[i + 1] <= age & num[i] < size_limit) {
      i <- i + 1
      if (verbose == TRUE) {
        start_time <- Sys.time()
      }
      ed <-
        edd_get_ed(num[i - 1], l_table, t[i], metric, offset, converter)
      lamu_real <- edd_update_lamu(ed, ed_max, params, model)
      if (verbose == TRUE) {
        message("Current rates:")
        print(lamu_real)
      }

      prob_real <- sum(lamu_real$newlas + lamu_real$newmus)

      prob_diff_la <- max(0, sum(lamu$newlas - lamu_real$newlas))
      prob_diff_mu <- max(0, sum(lamu$newmus - lamu_real$newmus))
      prob_fake <- sum(prob_diff_la + prob_diff_mu)

      if (verbose == TRUE) {
        message("Sum of probabilities of real events:")
        print(prob_real)
        message("Sum of probabilities of fake events:")
        print(prob_fake)
      }

      event_type <- sample(c("real", "fake"),
                           1,
                           prob = c(prob_real, prob_fake))

      if (verbose == TRUE) {
        message(paste0("Simulation step ", i, " started"))
        message(paste0("The next event happened is ", event_type))
      }

      if (event_type == "real") {
        event <-
          edd_sample_event(lamu_real$newlas, lamu_real$newmus, linlist)
        ran_lin <- c(linlist, linlist)[event]

        if (event <= length(linlist)) {
          num[i] <- num[i - 1] + 1
          new_lin <- new_lin + 1
          l_table <-
            rbind(l_table, c(t[i], ran_lin, sign(ran_lin) * new_lin, -1))
          linlist <- c(linlist, sign(ran_lin) * new_lin)

          if (verbose == TRUE) {
            message("Speciation event happened")
            message(paste0("Current species number is ", num[i]))
            message("Current species number is\n")
            print(l_table)
          }
        } else {
          num[i] <- num[i - 1] - 1
          l_table[abs(ran_lin), 4] <- t[i]
          w <- which(linlist == ran_lin)
          linlist <- linlist[-w]

          if (verbose == TRUE) {
            message("Extinction event happened")
            message(paste0("Current species number is ", num[i]))
            message("Current species number is\n")
            print(l_table)
          }
        }
      } else {
        num[i] <- num[i - 1]
      }

      if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
        t[i + 1] <- Inf

        if (verbose == TRUE) {
          message("One of the crown lineage goes extinct, simulation terminated")
        }
      } else {
        ed <- edd_get_ed(num[i], l_table, t[i], metric, offset, converter)
        ed_max <-
          edd_get_edmax(num[i], l_table, age, metric, offset, converter)
        params[1] <- num[i]
        lamu <- edd_update_lamu(ed, ed_max, params, model)

        if (verbose == TRUE) {
          message("Current rates:")
          print(lamu)
        }

        if (history == TRUE) {
          eds <- c(eds, list(ed))
          las <- c(las, list(lamu$newlas))
          mus <- c(mus, list(lamu$newmus))
          linlists <- c(linlists, list(linlist))
        }

        if (edd_sum_rates(lamu$newlas, lamu$newmus) == 0) {
          t[i + 1] <- Inf

          if (verbose == TRUE) {
            message("Sum of the rates was 0, simulation terminated")
          }
        } else {
          t[i + 1] <-
            t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))

          if (verbose == TRUE) {
            message(paste0("The next time interval is ", t[i + 1]))
          }
        }
      }

      if (verbose == TRUE) {
        end_time <- Sys.time()
        times[i] <- end_time - start_time
      }
    }

    if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
      done <- 0

      if (verbose == TRUE) {
        message("One of the crown lineage is extinct, simulation will start over")
      }
    } else if (num[i] >= size_limit) {
      done <- 0
      
      if (verbose == TRUE) {
        warning("Tree exceeded size limit, simulation will start over")
      }
    } else {
      done <- 1

      if (verbose == TRUE) {
        message("Simulation finished")
      }
    }
  }

  if (size_limit < 1e6) {
    tes <- l_to_phylo_ed(l_table, t[i] + 0.01, drop_extinct = TRUE)
    tas <- l_to_phylo_ed(l_table, t[i] + 0.01, drop_extinct = FALSE)
  } else {
    tes <- l_to_phylo_ed(l_table, age, drop_extinct = TRUE)
    tas <- l_to_phylo_ed(l_table, age, drop_extinct = FALSE)
  }

  ltt <-
    data.frame("time" = t[-length(t)],
               "num" = num)

  if (history == TRUE) {
    out <-
      list(
        tes = tes,
        tas = tas,
        l_table = l_table,
        ltt = ltt,
        eds = eds,
        las = las,
        mus = mus,
        linlists = linlists
      )
  } else {
    out <-
      list(
        tes = tes,
        tas = tas,
        l_table = l_table,
        ltt = ltt
      )
  }

  if (verbose == TRUE) {
    running_times <- data.frame("t" = times,
                                "n" = num)
    message("Results recorded")
    dir.create(file.path(getwd(), "/logs"), showWarnings = FALSE)
    utils::write.csv(running_times,
                     paste0(
                       getwd(),
                       "/logs/",
                       format(Sys.time(), "%Y-%m-%d_%H_%M_%S"),
                       ".csv"
                     ),
                     row.names = FALSE)
  }

  return(out)
}
