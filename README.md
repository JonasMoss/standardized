
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Please avoid the standardized alpha and the ordinal alpha<img src="man/figures/logo.png" align="right" width="140" height="70" />

[![Build
Status](https://travis-ci.com/JonasMoss/standardized.svg?branch=master)](https://travis-ci.org/JonasMoss/standardized)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/JonasMoss/standardized?branch=master&svg=true)](https://ci.appveyor.com/project/JonasMoss/standardized)

Repository for the PsyArXiv version of the paper “Please avoid the
standardized alpha and the ordinal alpha”.

## Overview

The new concrete ordinal reliabilities are implemented in a separate
package, [conogive](https://github.com/JonasMoss/conogive).

### Package

The `R` package of this paper includes functions for calculating
classical reliabilities, tables, et cetera.

### Paper

The paper itself is `standardized.tex`. Figures and Latex code generated
by `R` are in the `chunks` folder. The `R` folder contains functions
used in the paper or in the development of the paper. The folder `code`
contains the code used to generate the paper’s figures, simulations, and
so on.

## Installation

From inside `R`, use the following commands:

``` r
devtools::install_github("JonasMoss/standardized")
devtools::install_github("JonasMoss/conogive")
```

These packages needs to installed if you intend intend to run files in
the `code` directory.

## How to contribute

Read the paper\!
