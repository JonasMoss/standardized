
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Please avoid the standardized alpha <img src="man/figures/logo.png" align="right" width="140" height="70" />

[![Build
Status](https://travis-ci.org/JonasMoss/standardized.svg?branch=master)](https://travis-ci.org/JonasMoss/standardized)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/JonasMoss/standardized?branch=master&svg=true)](https://ci.appveyor.com/project/JonasMoss/standardized)
[![CircleCI build
status](https://circleci.com/gh/JonasMoss/standardized.svg?style=svg)](https://circleci.com/gh/JonasMoss/standardized)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/JonasMoss/standardized/branch/master/graph/badge.svg)](https://codecov.io/gh/JonasMoss/standardized?branch=master)

<!--[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)-->

<!--[![Project Status: Unsupported – The project has reached a stable, usable state but the author(s) have ceased all work on it. A new maintainer may be desired.](https://www.repostatus.org/badges/latest/unsupported.svg)](https://www.repostatus.org/#unsupported) -->

<!--[![DOI](https://zenodo.org/badge/120678148.svg)](https://zenodo.org/badge/latestdoi/120678148) -->

Joint repository for the unfinished paper “Please avoid the standardized
alpha” and its companion `R` package.

**N.B.** This repository will be split in two in the *near* future; one
package for the paper and one for the `R` package. The paper itself is
unfinished but should be preprint-ready soon, as by 15 March 2020.

## Overview

### Package

The package contains functions to estimate many different sorts of
reliabiltiy coefficients. The new ordinal reliability is in implemented
in the functions `ordinal_poly` (for estimation), `ordinal_omega`
(population value) and `ordinal_alpha` (the coefficient alpha variant.)

### Paper

The paper itself is in the `paper` folder. Figures and Latex code
generated by `R` are in the `chunks` folder. The `R` folder contains
functions used in the paper or in the development of the paper. The
folder `code` contains the code used to generate the paper’s figures,
simulations, and so on.

Some of the functions in this package could be of independent interest,
but I will not make it a `CRAN` package. Of course, feel free to grab
the code here for your own projects.

## Installation

From inside `R`, use one of the following commands:

``` r
devtools::install_github("JonasMoss/standardized")
```

The package needs to installed if you intend intend to run files in the
`code` directory.

## Usage Example

Run `run.R` from inside of `R` to repopulate the `chunks` directory.

## How to Contribute

I’m planning on being the lone author on this paper. Right now I don’t
welcome contributions.
