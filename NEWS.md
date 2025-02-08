# ribd 1.7.1

## New features
* `kappaIBD()` gains a logical argument `acrossComps`. Switching this off will compute the IBD coefficients only between individuals in the same connected component, which is often better for large pedigrees.

* Update and improve documentation of `multiPersonIBD()`.

* Several tweaks of the `plotly` version of IBD triangles.

## Other
* Now using Rhub v2 for testing.
* Increased pedtools dependency to 2.6.0, for `labels(..., unlist = TRUE)`.


# ribd 1.7.0

Th main focus of this version is an overhaul of the IBD triangle plots, and implementing such plots also in `ggplot2` and `plotly`.

## New features

* `ibdTriangle()` and `showInTriangle()` gain an argument `plotType`, with permitted values "base" (default), "ggplot2" and "plotly".

* `ibdTriangle()` now supports graphic parameters `las`, `mar` and `title`. (These work for base plots and also `ggplot2`.)

* `showInTriangle()` has a new argument `ped`, for adding an inset pedigree in the top right corner. 

* Improved handling of labels in `showInTriangle()`.

* New dataset `basicRelationships`, which is called by `ibdTriangle()`.

* In `kinship(x, ids)`, `ids` is no longer restricted to length 2. 

* Improve documentation of `kappaIBD()`.

## Bug fixes
* `kinship(x)` now catches duplicated ID labels across components when `x` is a pedlist.


# ribd 1.6.1

## New features

* `ibdTriangle()` gains argument `shortLines`, restricting kinship lines to the interior of the triangle.

* `ibdDraw()` now automatically calculates sensible plot margins.

* Add citation info.


## Bug fixes

* `kappaIBD(x, ids, simplify = TRUE)` now works as intended when `x` is a ped list and `ids` has length 2.

* Fixed printing of two-locus kinship patterns in `twoLocusIBD(..., verbose = TRUE)`.


# ribd 1.5.0

## Breaking changes

* In the output of `coeffTable()` several columns have been renamed.

## New features

* New function `twoLocusInbreeding()`.

* New function: `realisedIbdVariance()`.

* New function `ELR()`, implementing method of Egeland & Slooten (2016).

* All two-locus functions have been cleaned up and improved, and new examples have been added.

## Bug fixes

* Check that the input vector to `detailed2condensed()` has sum 1 (#14).

* Fixed edge-case bug in `twoLocusIdentity()` (#9).


# ribd 1.4.0

## Breaking changes

* The function `idcoefs()` was removed, since it relied on the no-longer-available package `identity`.

* The deprecated `kinshipX()`, `inbreedingX()` and `kappaIbdX()` were removed, and replaced with an argument `Xchrom` in `kinship()`, `inbreeding()` and `kappaIBD()`.

* The function `generalisedKinship()` has been replaced with the much more versatile `gKinship()` (see below). 

## New features

* New function `identityCoefs()` for computing condensed and detailed identity coefficients ("Jacquard coefficients"). Both autosomal and X-chromosomal versions are supported. This function supersedes `condensedIdentity()` and `condensedIdentityX()`, which will continue to exist, nonetheless. 

* `condensedIdentity()` gains arguments `simplify` and `self`, to match the new `identityCoefs()`.

* Computation of identity coefficients by MERLIN (via the "--extended" feature) is implemented in `identityCoefs()` with the option  `method = "merlin"`. Note that MERLIN rounds the output to 3 decimals, reducing its utility somewhat.

* New function `gKinship()` for computing generalised kinship coefficients. Several algorithms are implemented, supporting various flavours (random vs. deterministic; distinct vs. non-distinct groups; autosomal vs. X).

* New container class `gip` for generalised IBD patterns. Includes a print method.

* New function `coeffTable()` collecting various pedigree coefficients in a single table.

* New function `kin2deg()` computing the degree of relatedness, as used e.g. by the software KING.

## Bug fixes
* Fixed bug in `kinship(x, ids)` affecting pedigrees in nonstandard order.


# ribd 1.3.1

* `kinship()` now accepts a list of pedigrees as input.


# ribd 1.3.0

* The README has been rewritten and substantially expanded.

* `kappaIBD()` now accepts ped lists as input.

* In `inbreeding()` the argument `id` is renamed to `ids` and accepts vectors of length > 1.


# ribd 1.2.0

* New function `ibdDraw()` for illustrating IBD patterns in a pedigree. IBD alleles are represented as coloured dots or as letters.

* New function `ibdTriangle()`, which replaces `forrel::IBDtriangle()`.

* `constructPedigree()` now gives a textual description of the (usually double-half-cousin-like) pedigree it produces.


# ribd 1.1.0

* New function `constructPedigree()`, which constructs a pedigree yielding a prescribed set of IBD coefficients. This implements the algorithm described in https://doi.org/10.1007/s00285-020-01505-x.  

* `kinship()` gains a new argument `ids`, which is handy when you only want the kinship coefficient between two individuals (and not the whole kinship matrix)

* Similarly, `inbreeding()` gains the argument `id` for computing the inbreeding coefficient of a single individual.


# ribd 1.0.1

* This is a minor release, resolving a CRAN request to fix a documentation issue.

# ribd 1.0.0

* Initial CRAN release
