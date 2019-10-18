## Resubmission
This is a resubmission. In this version I have:

* As requested, removed code in examples that was commented out

* As requested, replaced calls to cat() with message(). One exception is in the function 
* A few other minor cleanups of examples

* Renamed functions 
jacquard() --> idcoefs() 
jacquard() --> idcoefs2()

I was asked by CRAN volunteer to explain the phrase "written by Abney (2009)" in the docs of these functions, since Abney is not a package author. 

These are wrapper functions. The first calls a function in the 'identity' package; the second attempts to call an external C program called "IdCoefs", if this happens to be installed on the users computer. This C program is written by Abney, which I've tried to make clear in the documentation. I'm the author of every line of code in the package.

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
