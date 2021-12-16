#' Interpret the files produced by my count_idx.pl script.
#'
#' This uses the index_table column in the metadata, reads this files associated
#' with it, and creates some data structures with the results.  These include a
#' modified version of the metadata, containing some summary information, a
#' table of the reads/tag observed, and a long table for ggplot.
#'
#' @param metadata Sample metadata.
#' @param id_column Column in the metadata containing the sample names.
#' @param index_column Column in the metadata containing the tag tables.
#' @param output Write out the matrix to this file (if provided).
#' @param cutoff Initial reads/tag filter.
#' @export
read_tags <- function(metadata, id_column = "sampleid", index_column = "index_table",
                     output = NULL, cutoff = 3) {
  long_samples <- data.frame()
  all_counts <- data.table::data.table()
  modified <- metadata
  modified[["pre_reads"]] <- 0
  modified[["post_reads"]] <- 0
  modified[["observed"]] <- 0
  modified[["passed_cutoff"]] <- 0
  modified[["rejected_cutoff"]] <- 0
  modified[["max_reads_per_tag"]] <- 0
  for (s in rownames(metadata)) {
    table <- metadata[s, index_column]
    sampleid <- metadata[s, id_column]
    element <- read.csv(file = table, sep = "\t", header = FALSE)
    colnames(element) <- c("idx", "number")
    modified[s, "pre_reads"] <- sum(element[["number"]])
    modified[s, "observed"] <- nrow(element)
    if (!is.null(cutoff)) {
      keep_idx <- element[["number"]] >= cutoff
      modified[s, "passed_cutoff"] <- sum(keep_idx)
      modified[s, "rejected_cutoff"] <- sum(!keep_idx)
      element <- element[keep_idx, ]
    }
    ## This tmp allocation is a little dumb, but for the moment
    ## I will use it to create a full matrix of barcodes.
    tmp <- data.table::as.data.table(element)
    colnames(tmp) <- c("idx", sampleid)
    if (ncol(all_counts) == 0) {
      all_counts <- tmp
    } else {
      all_counts <- merge(all_counts, tmp, by = "idx", all = TRUE)
    }
    element[["cfu"]] <- metadata[s, "cfu"]
    element[["sample"]] <- metadata[s, "condition"]
    element[["sampleid"]] <- sampleid
    modified[s, "post_reads"] <- sum(element[["number"]])
    modified[s, "max_reads_per_tag"] <- max(element[["number"]])
    long_samples <- rbind(long_samples, element)
  }

  retlist <- list(
    "modified_metadata" = modified,
    "long_samples" = long_samples,
    "sample_counts" = all_counts)
  if (!is.null(output)) {
    write.csv(x = all_counts, file = output)
  }
  return(retlist)
}
