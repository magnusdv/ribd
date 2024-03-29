% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realisedIbdVariance.R
\name{realisedIbdVariance}
\alias{realisedIbdVariance}
\title{Variance of realised relatedness coefficients}
\usage{
realisedIbdVariance(x, ids = leaves(x), coeff, L = 1)
}
\arguments{
\item{x}{A ped object.}

\item{ids}{A vector naming two members of \code{x}.}

\item{coeff}{A string naming a coefficient for which the variance is to be
computed. See Details for legal values.}

\item{L}{A positive number; the chromosome length in Morgan.}
}
\value{
A positive number.
}
\description{
Compute the variance of realised relatedness coefficients, by doubly
integrating the corresponding two-locus coefficients.
}
\details{
The double integral method was used by Guo to compute the variation in
proportion of the genome shared IBD (Guo 1995, see also Thompson 2013). The
method extends directly to other coefficients. The implementation here
supports Cotterman's kappa coefficients (of noninbred individuals), and
Jacquard's condensed identity coefficients.

This function is a bare-bones implementation of the double integral method,
based on \code{stats::integrate}, and can probably be optimised in various ways.

The \code{coeff} parameter must be either a character naming the coefficient to
compute, or a function. If a character, it must be one of the following
names:
\itemize{
\item "inb" (inbreeding coefficient)
\item "kinship", "phi" (synonyms for the kinship coefficient)
\item "k0", "k1", "k2" (kappa coefficients of noninbred individuals)
\item "D1", "D2", ... "D9" (condensed identity coefficients)
}
}
\examples{
###################################
### Box 1 of Hill & Weir (2011) ###
###################################

# Eq. 4b of Hill & Weir
phi = function(n, l) {
 1/(2*l^2) * (1/4)^n * sum(sapply(1:n, function(r)
   choose(n, r) * (2*r*l - 1 + exp(-2*r*l))/r^2))
}

# Chromosome of 1 Morgan
L = 1

### Full sibs ###

\dontrun{
x = nuclearPed(2)
realisedIbdVariance(x, ids = 3:4, coeff = "k2", L = L)

# Hill & Weir (Box 1)
16*phi(4,L) - 16*phi(3,L) + 8*phi(2,L) - 2*phi(1,L)
}

### Double first cousins ###

\dontrun{
dfc = doubleFirstCousins()

# Runtime ~1 min
realisedIbdVariance(dfc, coeff = "k0", L = L)
realisedIbdVariance(dfc, coeff = "k1", L = L)
realisedIbdVariance(dfc, coeff = "k2", L = L)

# Hill & Weir, Box 1
var_k2 = 64*phi(8,L) - 64*phi(7,L) + 40*phi(6,L) - 20*phi(5,L) +
  33/4*phi(4,L) - 5/2*phi(3,L) + 5/8*phi(2,L)-1/8*phi(1,L)
var_k1 = 4*var_k2
var_k0 = var_k2 + 2 * (4*phi(4,L) - 2*phi(3,L) + 3/4*phi(2,L) - 1/4*phi(1,L))

var_k0
var_k1
var_k2
}

}
\references{
\itemize{
\item Guo (1995) \emph{Proportion of genome shared identical by descent by relatives:
concept, computation, and applications}. Am J Hum Genet.
\item Hill & Weir (2011). \emph{Variation in actual relationship as a consequence of
Mendelian sampling and linkage}. Genet Res.
\item Thompson (2013). \emph{Identity by Descent: Variation in Meiosis, Across
Genomes, and in Populations}. Genetics.
}
}
