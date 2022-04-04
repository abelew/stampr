read_csv_tags <- function(input_csv, meta, id_column = "sampleid", cutoff = NULL,
                          sample_column = "replicate")  {
  input_df <- readr::read_csv(input_csv)
  colnames(input_df) <- tolower(gsub(pattern = "\\.txt", replacement = "",
                                     x = colnames(input_df)))
  colnames(input_df)[1] <- "idx"
  long_samples <- data.frame()
  modified <- meta
  modified[["csv_pre_reads"]] <- 0
  modified[["csv_post_reads"]] <- 0
  modified[["csv_observed"]] <- 0
  modified[["csv_passed_cutoff"]] <- 0
  modified[["csv_rejected_cutoff"]] <- 0
  modified[["csv_max_reads_per_tag"]] <- 0

  for (s in rownames(meta)) {
    ## I think the following line is redundant, the sample id is the rowname?
    sampleid <- meta[s, id_column]
    element <- input_df[, c("idx", sampleid)]
    colnames(element) <- c("idx", "number")
    modified[s, "csv_pre_reads"] <- sum(element[["number"]])
    modified[s, "csv_observed"] <- sum(element[["number"]] > 0)
    if (!is.null(cutoff)) {
      keep_idx <- element[["number"]] >= cutoff
      modified[s, "csv_passed_cutoff"] <- sum(keep_idx)
      modified[s, "csv_rejected_cutoff"] <- sum(!keep_idx)
      element <- element[keep_idx, ]
    }
    ## Now start appending to the long counts the observed reads/tag on a per-sample basis
    element[["cfu"]] <- meta[s, "cfu"]
    element[["sample"]] <- meta[s, sample_column]
    element[["sampleid"]] <- sampleid
    modified[s, "csv_post_reads"] <- sum(element[["number"]])
    modified[s, "csv_max_reads_per_tag"] <- max(element[["number"]], na.rm = TRUE)
    long_samples <- rbind(long_samples, element)
  }

  retlist <- list(
      "input_csv" = input_csv,
      "modified_metadata" = modified,
      "long_samples" = long_samples,
      "sample_counts" = input_df)
  return(retlist)
}
