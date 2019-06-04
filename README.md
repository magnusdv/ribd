<!-- README.md is generated from README.Rmd. Please edit that file -->

ribd <img src="man/figures/logo.png" align="right" height=140/>
===============================================================

Overview
--------

The goal of `ribd` is to compute various coefficients of relatedness
between pedigree members. It extends the `pedtools` package which
provides many useful utilities for pedigree construction and
manipulation.

The main functions in ribd are:

-   `kinship()`, `kinshipX()` : Computes kinship matrices
-   `inbreeding()`, `inbreedingX()` : Computes inbreeding coefficients
-   `kappaIbd()`, `kappaIbdX()` : Computes IBD coefficients
    (*κ*<sub>0</sub>, *κ*<sub>1</sub>, *κ*<sub>2</sub>) of non-inbred
    pedigree members
-   `condensedIdentity()`, `condensedIdentityX()` : Computes Jacquard’s
    condensed identity coefficients of any pedigree members

A novel feature of `ribd` is the ability to handle pedigrees with inbred
founders. More about this below!

The package also computes some more esoteric coefficients:

-   `generalisedKinship()` : Generalised kinship coefficients, as
    defined by Karigl (1981)
-   `twoLocusKinship()` : Two-locus kinship coefficients, as defined by
    Thompson (1988)
-   `twoLocusIBD()` : Two-locus IBD coefficients
-   `twoLocusGeneralisedKinship()` : Generalised two-locus kinship
    coefficients (*not exported*)

Installation
------------

`ribd` is available for download from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install ribd from github
devtools::install_github("magnusdv/ribd")
```

Getting started
---------------

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

Inbred founders
---------------

How would the above result would change if individual 1 was himself
inbred? The “normal” approach to answer this question would be to
*expand* the pedigree to include the complete family history. For
example, if he was the child of half siblings, we could perform this
expansion by merging `x` with a suitably labeled half-sib pedigree:

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
#> [1] 0.06640625 0.03515625 0.13281250 0.03125000 0.13281250 0.03125000
#> [7] 0.21875000 0.29687500 0.05468750
```

Although the above strategy worked nicely in this case, it quickly gets
awkward or impossible to model founder inbreeding by creating the
complete pedigree. For example, inbreeding coefficients close to zero
require enormous pedigrees! And even worse: What if individual 1 was
100% inbred? This cannot be modelled in this way, as it calls for an
infinite pedigree.

A much easier approach is to use the `founderInbreeding()` feature
offered by `pedtools`: We simply specify the inbreeding level of
individual 1 (in the original `x`) to be that of a child of half
siblings, i.e. 1/8.

``` r
founderInbreeding(x, ids = 1) = 1/8
```

When we now run `condensedIdentity()` on `x`, this inbreeding is taken
into account, giving the same answer as above.

``` r
condensedIdentity(x, ids = 5:6)
#> [1] 0.06640625 0.03515625 0.13281250 0.03125000 0.13281250 0.03125000
#> [7] 0.21875000 0.29687500 0.05468750
```

The condensed identity states
-----------------------------

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
combination, we re-use the autosomal states by using the following
simple rule: **Replace the single allele of any male, with a pair of
autozygous alleles**. In this way we get a one-to-one map from the X
states to the autosomal state.

For simplicity the output always contains 9 coefficients, but with NA’s
in the positions of undefined states (depending on the sex combination).
Hopefully this should all be clear from the following table:
<img src="man/figures/jacquardStatesX.png" align="left">
