#' @export edd_sim_nr
edd_sim_nr <- function(pars,
                    age,
                    model = "dsce2",
                    metric = "ed",
                    offset = "none",
                    history = TRUE,
                    verbose = FALSE,
                    converter = "cpp",
                    size_limit = 1e6) {
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
            message("Current L table is\n")
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
            message("Current L table is\n")
            print(l_table)
          }
        }
      } else {
        num[i] <- num[i - 1]

        if (verbose == TRUE) {
          message("Fake event happened, nothing changed")
          message(paste0("Current species number is ", num[i]))
          message("Current L table is\n")
          print(l_table)
        }
      }

      if (num[i] < 2) {
        t[i + 1] <- Inf

        if (verbose == TRUE) {
          message("Tree had tips less than 2, simulation will terminated")
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

    if (num[i] < 2) {
      done <- 0

      if (verbose == TRUE) {
        message("Tree had tips less than 2, simulation will restart")
      }
    } else {
      done <- 1

      if (verbose == TRUE) {
        message("Simulation successfully finished")
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
    message("Simulation time steps recorded")
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