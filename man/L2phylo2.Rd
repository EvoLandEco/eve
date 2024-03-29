% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edd_sim_utilities.R
\name{L2phylo2}
\alias{L2phylo2}
\title{Function to convert a table with speciation and extinction events to a
phylogeny}
\usage{
L2phylo2(L, t, dropextinct = FALSE)
}
\arguments{
\item{L}{Matrix of events as produced by dd_sim: \cr \cr - the first column
is the time at which a species is born in Mya\cr - the second column is the
label of the parent of the species; positive and negative values indicate
whether the species belongs to the left or right crown lineage \cr - the
third column is the label of the daughter species itself; positive and
negative values indicate whether the species belongs to the left or right
crown lineage \cr - the fourth column is the time of extinction of the
species; if the fourth element equals -1, then the species is still extant.}

\item{t}{Simulation time}

\item{dropextinct}{Sets whether the phylogeny should drop species that are
extinct at the present}
}
\value{
\item{ phy }{ A phylogeny of the phylo type }
}
\description{
Converting a table with speciation and extinction events to a phylogeny
(reversed time scale)
}
\examples{
# do not use this function, use L2phylo() instead

}
\references{
\itemize{
\item Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
180: E75-E89, doi: 10.1086/667574
}
}
\author{
Rampal S. Etienne; Tianjian Qin
}
\keyword{models}
