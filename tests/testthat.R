#!/usr/bin/env Rscript
library(testthat)
library(stampr)

test_package("stampr", reporter="summary", stop_on_failure=FALSE)
