% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edd_wrapper.r
\name{edd_sim_rep}
\alias{edd_sim_rep}
\title{edd_sim_rep}
\usage{
edd_sim_rep(combo = NULL, history = FALSE, verbose = FALSE, nrep = 5)
}
\arguments{
\item{combo}{A list containing parameter sets to be used to start a simulation}

\item{history}{Logical, decides whether to record historical states for all
of the lineages, including the transitions of the rates and the evolutionary
relationship metrics}

\item{verbose}{Logical, decides whether to print simulation details of each}

\item{nrep}{number of replication
time step during simulation, for debugging purposes}
}
\description{
edd_sim_rep
}
\author{
Tianjian Qin
}
