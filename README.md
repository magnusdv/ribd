
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ribd <img src="man/figures/logo.png" align="right" height=140/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ribd)](https://CRAN.R-project.org/package=ribd)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ribd?color=yellow)](https://cran.r-project.org/package=ribd)
[![](https://cranlogs.r-pkg.org/badges/last-month/ribd?color=yellow)](https://cran.r-project.org/package=ribd)
<!-- badges: end -->

## Overview

The goal of `ribd` is to compute various coefficients of relatedness and
identity-by-descent (IBD) between pedigree members. It extends the
`pedtools` package which provides useful utilities for pedigree
construction and manipulation.

The main functions in `ribd` are the following:

  - `kinship()`, `kinshipX()` : Kinship coefficients
  - `inbreeding()`, `inbreedingX()` : Inbreeding coefficients
  - `kappaIBD()`, `kappaIBDX()` : IBD coefficients
    \((\kappa_0, \kappa_1, \kappa_2)\) (noninbred individuals only)
  - `condensedIdentity()`, `condensedIdentityX()` : Jacquard’s condensed
    identity coefficients

A unique feature of `ribd` is the ability to handle pedigrees with
inbred founders in all of the above calculations. More about this
below\!

The package also computes a variety of lesser-known pedigree
coefficients:

  - `generalisedKinship()` : Generalised kinship coefficients, as
    defined by Weeks and Lange (1988)
  - `multiPersonIBD()` : Multi-person IBD coefficients (noninbred
    individuals only)
  - `twoLocusKinship()` : Two-locus kinship coefficients, as defined by
    Thompson (1988)
  - `twoLocusIBD()` : Two-locus IBD coefficients (noninbred pair of
    individuals)
  - `twoLocusIdentity()` : Two-locus condensed identity coefficients
    (any pair of individuals)
  - `twoLocusGeneralisedKinship()` : Generalised two-locus kinship
    coefficients (*not exported*)

## Installation

To get the current official version of `ribd`, install from CRAN as
follows:

``` r
install.packages("ribd")
```

Alternatively, you can obtain the latest development version from
GitHub:

``` r
# install.packages("devtools") # install devtools if needed
devtools::install_github("magnusdv/ribd")
```

## Getting started

``` r
library(ribd)
```

To illustrate the use of `ribd` we compute the condensed identity
coefficients after one generation of full sib mating. This is a suitable
example because the answer is well known, and it is one of the simplest
in which all 9 coefficients are non-zero.

We create the pedigree with `pedtools` as follows:

``` r
x = fullSibMating(1)
plot(x)
```

<img src="man/figures/README-sibs-1.png" style="display: block; margin: auto;" />

The identity coefficients of the children are computed with
`condensedIdentity()`

``` r
condensedIdentity(x, ids = 5:6)
#> [1] 0.06250 0.03125 0.12500 0.03125 0.12500 0.03125 0.21875 0.31250 0.06250
```

## Inbred founders

How would the above result would change if individual 1 was himself
inbred, say, as a child of half siblings? A possible, but cumbersome,
approach to answer this question would be to *expand* the pedigree to
include the complete family history. Here is one way to do this, using
`pedtools::mergePed()` to merge `x` with a suitably labelled half-sib
pedigree:

``` r
y = halfSibPed(sex1 = 1, sex2 = 2)
y = addChildren(y, father = 4, mother = 5, nch = 1)
y = relabel(y, c(101:105, 1)) # prepare merge by relabeling
z = mergePed(y, x)
```

``` r
plot(z)
```

<img src="man/figures/README-sibs-extended-1.png" style="display: block; margin: auto;" />

Now that we have the complete pedigree we could answer the original
question by running `condensedIdentity()` on `z`.

``` r
condensedIdentity(z, ids = 5:6)
#> [1] 0.06640625 0.03515625 0.13281250 0.03125000 0.13281250 0.03125000 0.21875000
#> [8] 0.29687500 0.05468750
```

Although the above strategy worked nicely in this case, it quickly gets
awkward or impossible to model founder inbreeding by creating the
complete pedigree. For example, inbreeding coefficients close to zero
require enormous pedigrees\! And even worse: What if individual 1 was
100% inbred? This cannot be modelled in this way, as it calls for an
infinite pedigree.

A much easier approach is to use the `founderInbreeding()` feature
offered by `pedtools`: We simply specify the inbreeding level of
individual 1 (in the original `x`) to be that of a child of half
siblings, i.e. \(1/8\).

``` r
founderInbreeding(x, ids = 1) = 1/8
```

When we now run `condensedIdentity()` on `x`, this inbreeding is taken
into account, giving the same answer as for `z` above.

``` r
condensedIdentity(x, ids = 5:6)
#> [1] 0.06640625 0.03515625 0.13281250 0.03125000 0.13281250 0.03125000 0.21875000
#> [8] 0.29687500 0.05468750
```

## The pairwise condensed identity states

The following figure shows the 9 *condensed identity states* of two
individuals *a* and *b*. Each state shows the pattern of IBD between the
4 homologous alleles at an autosomal locus. The states are shown in the
ordering used by Jacquard and most subsequent authors. This is also the
order of the coefficients output by `condensedIdentity()`.
<img src="man/figures/jacquardStates.png" align="left">

### Identity states on X

The X chromosomal version of `condensedIdentity()` is called
`condensedIdentityX()`. What this function computes requires some
explanation, which we offer here.

The point of `condensedIdentityX()` is to compute the coefficients
(i.e., the expected proportions) of the pairwise identity states for a
locus on the X chromosome. What these states are depends on the sexes of
the involved individuals: either female-female, female-male, male-female
or male-male. In some sense the first case is the easiest: When both are
female the states are the same as in the autosomal case.

Males, being hemizygous, have only 1 allele of a locus on X. Hence when
males are involved the total number of alleles is less than 4, rendering
the autosomal states pictured above meaningless. However, to avoid
drawing (and learning the ordering of) new states for each sex
combination, we can re-use the autosomal state pictograms by invoking
the following simple rule: **Replace the single allele of any male, with
a pair of autozygous alleles**. This gives a one-to-one map from the X
states to the autosomal states.

For simplicity the output always contains 9 coefficients, but with NA’s
in the positions of undefined states (depending on the sex combination).
Hopefully this should all be clear from the following table:
<img src="man/figures/jacquardStatesX.png" align="left">
