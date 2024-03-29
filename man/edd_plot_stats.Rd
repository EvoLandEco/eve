% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edd_plot.r
\name{edd_plot_stats}
\alias{edd_plot_stats}
\title{Generating boxplots of tree statistics}
\usage{
edd_plot_stats(
  raw_data = NULL,
  method = "treestats",
  save_plot = FALSE,
  path = NULL
)
}
\arguments{
\item{raw_data}{a list of results generated by edd simulation function}

\item{method}{Specify which package to be used to calculate calculate statistics}

\item{save_plot}{Logical, whether save to file or return a ggplot object}

\item{path}{The path to save the plots}
}
\value{
an ggplot object
}
\description{
Function to generate boxplots showing the tree statistics e.g. Colless, J-One, B1, Gamma etc.
}
\author{
Tianjian Qin
}
\keyword{phylogenetics}
