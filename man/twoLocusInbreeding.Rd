% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoLocusInbreeding.R
\name{twoLocusInbreeding}
\alias{twoLocusInbreeding}
\title{Two-locus inbreeding}
\usage{
twoLocusInbreeding(x, id, rho, verbose = FALSE, debug = FALSE)
}
\arguments{
\item{x}{A pedigree in the form of a \code{\link[pedtools:ped]{pedtools::ped}} object.}

\item{id}{The ID label of a pedigree member.}

\item{rho}{A numeric vector of recombination rates; all entries must be in
the interval \eqn{[0, 0.5]}.}

\item{verbose}{A logical.}

\item{debug}{A logical. If TRUE, detailed messages are printed during the
recursion process.}
}
\description{
Computes the two-locus inbreeding coefficient of a pedigree member, for a
given recombination rate.
}
\details{
Let A be a pedigree member, and L1, L2 two autosomal loci with recombination
rate \eqn{\rho}. The two-locus inbreeding coefficient \eqn{f_{11}(\rho)} is
defined as the probability that A is autozygous at both L1 and L2
simultaneously.

As in the one-locus case, the two-locus inbreeding coefficient of A equals
the two-locus kinship coefficient of the parents.
}
\examples{

###################################################
# Reproducing an example of Weir & Cockerham (1969)
###################################################

# Pedigree
x = nuclearPed(2, sex = 1:2) |>
  addDaughter(3:4) |>
  addSon(c(3,5)) |>
  addDaughter(5:6) |>
  relabel(new = strsplit("GHDECBA","")[[1]])

plot(x)

# The two-locus inbreeding of A
twoLocusPlot(list(ped = x, ids = "A"), coeff = "inb")

# W&C formula (expressed by linkage parameter a = 1-2*rho)
rho = seq(0, 0.5, length = 11)
a = 1 - 2*rho
WC = (128 + 10*a + 36*a^2 + 47*a^3 + 20*a^4 + 10*a^5 + 4*a^6 + a^7)/512

points(rho, WC, col = 2)

}
\references{
Weir & Cockerham (1969). \emph{Pedigree mating with two linked loci}.
Genetics, 61:923-940.
}
\seealso{
\code{\link[=twoLocusKinship]{twoLocusKinship()}}, \code{\link[=twoLocusIBD]{twoLocusIBD()}}, \code{\link[=twoLocusIdentity]{twoLocusIdentity()}}
}
