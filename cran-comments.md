## Resubmission
This is a resubmission. In this version I have:

* As requested, replaced all T/F with TRUE/FALSE


## Resubmission
This is a resubmission. In this version I have:

* As requested, removed code in examples that was commented out.

* Some additional very minor cleanups of examples

* As requested, replaced cat()/print() with message() several places

* Renamed functions jacquard() --> idcoefs() and jacquard2() --> idcoefs2()

I was asked by CRAN volunteer to explain the phrase "written by Abney (2009)" in the docs of these functions, since Abney is not a package author. 

These are wrapper functions. The first calls a function in the 'identity' package; the second attempts to call an external C program called "IdCoefs", if this happens to be installed on the users computer. This C program is written by Abney, and is not shipped with the submitted package. I've tried to make this clear in the documentation. I'm the author of every line of code in the submitted package.


## Resubmission
This is a resubmission. In this version I have:

* Used "Title Case" for the package title.

## Test environments
* local Windows 10 install, R 3.6.1
* rhub::check_for_cran()
* devtools::check_win_devel()

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
