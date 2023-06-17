edd_rtree <- function(pars, age, model, metric, offset, size = NULL, n = 1,
                      history = FALSE, verbose = FALSE) {
  if (!(n %% 1 == 0)) {
    stop("n must be an integer")
  }

  if (!is.null(size)) {
    if (!(size %% 1 == 0)) {
      stop("size must be an integer")
    }
    if (size <= 2) {
      stop("size must be greater than 2")
    }
  }

  if (is.null(size)) {
    message("Searching for ", n, " trees, no size limit")
  } else {
    message("Searching for ", n, " trees with ", size, " tips")
  }

  edd_message_info(pars = pars, age = age, model = model, metric = metric, offset = offset)

  progressr::handlers(list(
    progressr::handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 60,
      complete = "+"
    )
  ))

  progress_search <-
    progressr::progressor(steps = n)

  trees <- list()
  i <- 0

  while(length(trees) < n) {
    rs <- edd_sim(pars = pars,
                  age = age,
                  model = model,
                  metric = metric,
                  offset = offset,
                  history = history,
                  verbose = verbose)
    if (is.null(size)) {
      i <- i + 1
      trees[[i]] <- rs$tes
      progress_search()
    } else {
      if (rs$tes$Nnode == size - 1) {
        i <- i + 1
        trees[[i]] <- rs$tes
        progress_search()
      }
    }
  }

  return(trees)
}


ddd_rtree <- function(pars, age, ddmodel, size = NULL, n = 1) {
  if (!(n %% 1 == 0)) {
    stop("n must be an integer")
  }

  if (!is.null(size)) {
    if (!(size %% 1 == 0)) {
      stop("size must be an integer")
    }
    if (size <= 2) {
      stop("size must be greater than 2")
    }
  }

  if (is.null(size)) {
    message("Searching for ", n, " trees, no size limit")
  } else {
    message("Searching for ", n, " trees with ", size, " tips")
  }

  progressr::handlers(list(
    progressr::handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 60,
      complete = "+"
    )
  ))

  progress_search <-
    progressr::progressor(steps = n)

  trees <- list()
  i <- 0

  while(length(trees) < n) {
    rs <- DDD::dd_sim(pars = pars,
                      age = age,
                      ddmodel = ddmodel)
    if (is.null(size)) {
      i <- i + 1
      trees[[i]] <- rs$tes
      progress_search()
    } else {
      if (rs$tes$Nnode == size - 1) {
        i <- i + 1
        trees[[i]] <- rs$tes
        progress_search()
      }
    }
  }

  return(trees)
}


#' @title Function to search for trees with a given number of tips and a given parameter set
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
#' @param size Specifies the number of tips in the tree to be searched for
#' @param n Specifies the number of trees to be searched for
#' @param history Logical, indicating whether to record the historical states
#' (of the rates and ED/PD values)
#' @param verbose Logical, for debugging purpose, indicating whether to print
#' simulation info at each step in the console, and save running time to a file
#' @param converter Either "cpp" or "r", choose which version of L2phylo to use.
#' @return a list of trees
#' @author Tianjian Qin
#' @export search_tree
search_tree <- function(pars, age, model, metric, offset, size = NULL, n = 1,
                        history = FALSE, verbose = FALSE) {
  progressr::with_progress({
    edd_rtree(pars = pars, age = age, model = model, metric = metric, offset = offset, size = size, n = n,
              history = history, verbose = verbose)
  })
}