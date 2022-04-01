
#' Plot the calibration curve of the tags observed vs. CFU
#'
#' Currently this is just a simple ggplot scatterplot.  I want to make it more
#' configurable and fun.
#'
#' @param tags Tag data.
#' @param x_column Metadata column for the x axis.
#' @param y_column Metadata column for the y axis.
#' @param color_column Metadata column for colors.
#' @param transform_x Perform a transformation on the x axis?
#' @param transform_y Perform a transformation on the y axis?
#' @export
exclude_samples <- function(tags, max_dilution = NULL, max_cfu = NULL, type = "calibration",
                            sampleids = NULL, replicate = NULL)  {
  meta <- tags[["modified_metadata"]]
  query_idx <- meta[["type"]] == type
  test_samples <- meta[query_idx, ]
  tag_data <- tags[["sample_counts"]]

  if (!is.null(max_dilution)) {
    message("Before dilution filter, there were: ", nrow(test_samples), " ", type, " samples.")
    keep_idx <- test_samples[["dilution"]] <= max_dilution
    test_samples <- test_samples[keep_idx, ]
    kept_samplenames <- rownames(test_samples)
    tag_data <- tag_data[, kept_samplenames]
    message("After dilution filter, there are: ", nrow(test_samples), " ", type, "samples.")
  }

  if (!is.null(max_cfu)) {
    message("Before cfu filter, there were: ", nrow(test_samples), " ", type, " samples.")
    keep_idx <- test_samples[["cfu"]] <= max_cfu
    test_samples <- test_samples[keep_idx, ]
    kept_samplenames <- rownames(test_samples)
    tag_data <- tag_data[, kept_samplenames]
    message("After cfu filter, there are: ", nrow(test_samples), " ", type, "samples.")
  }

  if (!is.null(sampleids)) {
    message("Before sampleid filter, there were: ", nrow(test_samples), " ", type, " samples.")
    keep_idx <- ! rownames(test_samples) %in% sampleids
    test_samples <- test_samples[keep_idx, ]
    kept_samplenames <- rownames(test_samples)
    tag_data <- tag_data[, kept_samplenames]
    message("After sampleid filter, there are: ", nrow(test_samples), " ", type, "samples.")
  }

  if (!is.null(replicate)) {
    message("Before replicate filter, there were: ", nrow(test_samples), " ", type, " samples.")
    keep_idx <- test_samples[["replicate"]] != replicate
    test_samples <- test_samples[keep_idx, ]
    kept_samplenames <- rownames(test_samples)
    tag_data <- tag_data[, kept_samplenames]
    message("After replicate filter, there are: ", nrow(test_samples), " ", type, "samples.")
  }

  retlist <- list(
      "previous_metadata" = tags[["previous_metadata"]],
      "modified_metadata" = test_samples,
      "Nb_column" = tags[["Nb_column"]],
      "sample_counts" = tag_data)
  return(retlist)
}
