
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eve <a href='https://github.com/EvoLandEco/eve/'><img src='man/eve-logos/eve-logos_transparent.png' align="right" height="139" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/furrr)](https://cran.r-project.org/package=furrr)
[![R build
status](https://github.com/DavisVaughan/furrr/workflows/R-CMD-check/badge.svg)](https://github.com/DavisVaughan/furrr/actions)
[![Codecov test
coverage](https://codecov.io/gh/DavisVaughan/furrr/branch/master/graph/badge.svg)](https://codecov.io/gh/DavisVaughan/furrr?branch=master)
<!-- badges: end -->

## Overview

The package `eve` is an evolution emulator which provides pipelines to
do phylogenetic-diversity-dependent simulation, analyse outputs and
generate publication-ready plots and tables conveniently. It serves as a
companion package to [DDD](https://github.com/rsetienne/DDD) and
[DAISIE](https://github.com/rsetienne/DAISIE) now, and maybe expanded to
accommodate more models and functions.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("rsetienne/DDD@tianjian_Rampal")
remotes::install_github("EvoLandEco/eve")
```

## Example

`eve` now only supports EDD simulation from the `DDD` package. The first
thing to do is to generate a combo consists of all the parameter sets we
need:

``` r
library(eve)

# make a combo of parameter sets
combo <- edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0.1, 0.2),
  beta_n = c(-0.001, 0),
  beta_phi = c(-0.001, 0.001),
  gamma_n = c(-0.001, 0.001),
  gamma_phi = c(-0.001, 0.001),
  age = 5,
  model = "dsde2",
  metric = c("ed"),
  offset = "none"
)

# have a look at the combo
combo[1:3]
#> $`1`
#>   age model metric offset                                         pars
#> 1   5 dsde2     ed   none 0.500, 0.100, -0.001, -0.001, -0.001, -0.001
#> 
#> $`2`
#>   age model metric offset                                         pars
#> 2   5 dsde2     ed   none 0.800, 0.100, -0.001, -0.001, -0.001, -0.001
#> 
#> $`3`
#>   age model metric offset                                         pars
#> 3   5 dsde2     ed   none 0.500, 0.200, -0.001, -0.001, -0.001, -0.001
```

A large combo may take too long to finish, as an example a minimal combo
will do:

``` r
# make a minimal combo
combo <- edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0.1, 0.2),
  beta_n = -0.001,
  beta_phi = -0.001,
  gamma_n = 0.001,
  gamma_phi = 0.001,
  age = 5,
  model = "dsde2",
  metric = "ed",
  offset = "none"
)
```

`eve` supports parallel simulation, the following example shows how to
do a 8-session parallel EDD simulation given the combo, with 3
replications for each parameter set . The result will be saved to
`/result/example`. If `name` is not specified, a folder will be created
according to time and date at the moment.

``` r
# result will be save as .RData file
edd_go(
  combo = combo,
  nrep = 3,
  name = "example",
  strategy = future::multisession,
  workers = 8
)
#> C:/path/result/example already exists
#> Folder created
#> Running multisession simulation with 8 workers
#> Size of parameter space is: 4
#> Number of replications for each parameter set is: 3
#> 
#> Saving result to C:/path/result/example/example.RData
#> Result saved
```

`edd_go` will not return any object unless you set `name = "no_save"`.

``` r
# set name = "no_save" if you want to assign the output to a variable
output <- edd_go(
  combo = combo,
  nrep = 3,
  name = "no_save",
  strategy = future::multisession,
  workers = 8
)

# examine the result
output[[1]]
```

Sequential simulation is also possible:

``` r
edd <- edd_go(
    combo = combo,
    nrep = 3,
    name = "example2",
    strategy = future::sequential,
    workers = 8
  )
#> C:/path/result/example2 already exists
#> Folder created
#> Running sequential simulation
#> Size of parameter space is: 4
#> Number of replications for each parameter set is: 3
#> 
#> Saving result to C:/path/result/example2/example2.RData
#> Result saved
```
