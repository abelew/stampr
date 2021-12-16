#' Try out some heuristics for filtering the tag data.
#'
#' Given our focus on RNAseq, it is not a stretch to guess that we would assume
#' some of those methods would be useful.  Here is a cpm implementation of some
#' simple filtering.
#'
#' @param tag_list Result from one form of read_tags().  This contains the
#'  metadata and the matrix of reads/tag.
#' @param multiplier Arbitrary multiplier to make the filter more stringent.
#' @export
filter_tags <- function(tag_list, multiplier = 1) {
  modified <- tag_list[["modified_metadata"]]
  ## I have lots of ideas about filtering, here is the simplest
  all_mtrx <- as.data.frame(tag_list[["sample_counts"]])
  rownames(all_mtrx) <- all_mtrx[["idx"]]
  all_mtrx[["idx"]] <- NULL
  na_idx <- is.na(all_mtrx)
  all_mtrx[na_idx] <- 0
  all_mtrx <- edgeR::cpm(all_mtrx)
  ## Here is an aggressive possibility for filtering
  null_rows <- rowSums(all_mtrx) < (ncol(all_mtrx) * multiplier)
  aggressive <- all_mtrx[!null_rows, ]

  modified[["passed_aggressive_cpm"]] <- 0
  for (col in colnames(all_mtrx)) {
    agg_idx <- aggressive[, col] > 0
    modified[col, "passed_aggressive_cpm"] <- sum(agg_idx)
  }
  retlist <- list(
    "filtered" = aggressive,
    "modified_metadata" = modified)
  return(retlist)
}

#' Keep only the top-n most abundant tags
#'
#' In Dr. Lee's experiment, they explicitly know the number of input tags.  So
#' we can explicitly keep only those n most observed tags.
#'
#' @param tag_list Existing tag data
#' @param replicate_column Tags should be shared across replicates.
#' @param topn How many tags to keep.
#' @export
filter_topn_tags <- function(tag_list, replicate_column = "replicate", topn = 600) {
  modified <- tag_list[["modified_metadata"]]
  ## I have lots of ideas about filtering, here is the simplest
  all_mtrx <- as.data.frame(tag_list[["sample_counts"]])
  rownames(all_mtrx) <- all_mtrx[["idx"]]
  all_mtrx[["idx"]] <- NULL
  original_order <- colnames(all_mtrx)
  final_mtrx <- data.frame()
  for (replicate in levels(as.factor(modified[[replicate_column]]))) {
    rep_idx <- modified[[replicate_column]] == replicate
    rep_mtrx <- all_mtrx[, rep_idx]
    rep_sum <- rowSums(rep_mtrx, na.rm = TRUE)
    rep_order <- order(rep_sum, decreasing = TRUE)
    rep_ordered <- rep_mtrx[rep_order, ]
    desired_tags <- head(rownames(rep_ordered), n = topn)
    desired_mtrx <- rep_mtrx[desired_tags, ]
    if (ncol(final_mtrx) == 0) {
      final_mtrx <- desired_mtrx
    } else {
      final_mtrx <- merge(final_mtrx, desired_mtrx, by = "row.names", all = TRUE)
      rownames(final_mtrx) <- final_mtrx[["Row.names"]]
      final_mtrx[["Row.names"]] <- NULL
    }
  }
  ## Make certain that we did not reorder the matrix.
  final_mtrx <- final_mtrx[, original_order]
  modified <- modified[original_order, ]
  colname <- glue::glue("topn{topn}_totalreads")
  modified[[colname]] <- colSums(final_mtrx, na.rm = TRUE)
  colname <- glue::glue("topn{topn}_meanreads")
  modified[[colname]] <- colMeans(final_mtrx, na.rm = TRUE)
  colname <- glue::glue("topn{topn}_medianreads")
  modified[[colname]] <- matrixStats::colMedians(as.matrix(final_mtrx), na.rm = TRUE)
  tag_list[["modified_metadata"]] <- modified
  tag_list[["sample_counts"]] <- final_mtrx
  tag_list[["long_samples"]] <- NULL
  return(tag_list)
}
