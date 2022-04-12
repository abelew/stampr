## The following was taken from ggplot2's ggplot2.r
## I presume it is a blanket importer cue for roxygen2 to add
## import statements to the NAMESPACE file so that when ggplot2 is used
## it will ensure that these libraries are available.
## I checked the roxygen documentation and it appears that
## imports are saved as the exclusive set, as a result repeating these
## at each function declaration serves to make explicit what each function
## requires while not (I think) adding excessive cruft to the NAMESPACE

#' Pipe operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' dopar
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %dopar%
#' @rdname dopar
#' @keywords internal
#' @export
#' @importFrom foreach %dopar%
NULL

#' stampr: Some functions for analyzing sequencing data intended to quantify RT error rates.
#'
#' To see examples of this in action, check out the vignettes:
#' \code{browseVignettes(package = 'stampr')}
#'
#' @docType package
#' @name stampr
#' @importFrom glue glue glue_data
#' @importFrom foreach foreach
#' @import methods data.table dplyr ggplot2
#' @import utils
#' @import purrr
NULL
