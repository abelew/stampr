#' Try out some methods to go from tags/Nb back to CFU.
#'
#' This function has not been completed, mostly because I am not sure if it will
#' ever be used, but also because I am not sure of the appropriate model.  I
#' have never had an opportunity to use modelling to estimate data, so this is
#' mostly a place for me to play around with something that I should know but
#' don't.
#'
#' @param meta Metadata.
#' @param from Factor from which to extrapolate.
#' @param to Factor to which to return.
#' @param provided Dataframe of incomplete data.
#' @export
predict_cfu <- function(meta, from = "Nb", to = "cfu", provided = NULL) {
  calibration_idx <- meta[["type"]] == "calibration"
  cali <- meta[calibration_idx, ]
  cali[["logfrom"]] <- log10(cali[[from]])
  cali[["logto"]] <- log10(cali[[to]])
  lm_result <- glm(formula = logto ~ logfrom, data = cali)

  test_data <- data.frame()
  provided <- c(1.7, 2.4, 3.4, 4.3)
  if (!is.null(provided)) {
    test_data <- data.frame(row.names = 1:length(provided))
    test_data[["logfrom"]] = provided
  }
  prediction <- predict.glm(object = lm_result, newdata = test_data)
}
