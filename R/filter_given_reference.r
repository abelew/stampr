#' Filter samples using the tags observed in the reference samples.
#'
#' Take the result from read_tags() and use it to find which tags are in the
#' non-reference samples.  One caveat about this process, it splits the data by
#' replicate because tags are not global.  Thus we will need to consider what is
#' the best way to handle the filtered data.  The simplest is to just merge the
#' pieces back together.
#'
#' @param index_list Result from count_tags
#' @param reference_cutoff Minimum number of reads in the reference samples to
#'  be considered 'real'
#' @param remove_nas Replace NA with 0?
#' @export
filter_given_reference <- function(tag_list, reference_cutoff = 1, remove_nas = FALSE) {
  modified <- tag_list[["modified_metadata"]]
  final_mtrx <- as.data.frame(tag_list[["sample_counts"]])
  rownames(final_mtrx) <- final_mtrx[["idx"]]
  final_mtrx[["idx"]] <- NULL
  original_order <- colnames(final_mtrx)
  filtered_list <- filter_given_reference_list(tag_list, reference_cutoff)
  filtered_dt <- data.table::as.data.table(filtered_list[[1]], keep.rownames = TRUE)
  for (r in 2:length(filtered_list)) {
    tbl <- data.table::as.data.table(filtered_list[[r]], keep.rownames = TRUE)
    filtered_dt <- merge(filtered_dt, tbl, by = "rn", all = TRUE)
  }
  filtered_df <- as.data.frame(filtered_dt)
  rownames(filtered_df) <- filtered_dt[["rn"]]
  filtered_df[["rn"]] <- NULL
  if (isTRUE(remove_nas)) {
    na_idx <- is.na(filtered_df)
    filtered_df[na_idx] <- 0
  }
  final_mtrx <- filtered_df
  print(colnames(final_mtrx))
  print(original_order)
  final_mtrx <- final_mtrx[, original_order]
  modified <- modified[original_order, ]
  colname <- glue::glue("reffilt_totalreads")
  modified[[colname]] <- colSums(final_mtrx, na.rm = TRUE)
  colname <- glue::glue("reffilt_meanreads")
  modified[[colname]] <- colMeans(final_mtrx, na.rm = TRUE)
  colname <- glue::glue("reffilt_medianreads")
  modified[[colname]] <- matrixStats::colMedians(as.matrix(final_mtrx), na.rm = TRUE)
  tag_list[["modified_metadata"]] <- modified
  tag_list[["sample_counts"]] <- final_mtrx
  tag_list[["long_samples"]] <- NULL

  return(tag_list)
}

#' This does the real work for filter_given_reference().
#'
#' Create a list of matrices containing only the tags observed in each
#' replicates' reference sample.
#'
#' @param index_list Result from count_tags
#' @param reference_cutoff Minimum number of reads in the reference samples to
#'  be considered 'real'
filter_given_reference_list <- function(index_list, reference_cutoff = 1) {
  meta <- index_list[["modified_metadata"]]
  reads_per_tag <- index_list[["sample_counts"]]
  tags_per_replicate <- gather_tags_per_replicate(meta, reads_per_tag,
                                                  reference_cutoff = reference_cutoff)
  reference_sample_idx <- meta[["type"]] == "reference"
  reference_samples <- meta[reference_sample_idx, ]
  reference_samplenames <- rownames(reference_samples)
  reference_df <- as.data.frame(reads_per_tag)[, reference_samplenames]
  rownames(reference_df) <- reads_per_tag[["idx"]]

  experimental_sample_idx <- !reference_sample_idx
  experimental_samples <- meta[experimental_sample_idx, ]
  experimental_samplenames <- rownames(experimental_samples)
  experimental_df <- as.data.frame(reads_per_tag)[, experimental_samplenames]
  rownames(experimental_df) <- reads_per_tag[["idx"]]

  replicate_list <- list()
  message("Starting from ", nrow(experimental_df), " tags.")
  for (rep in names(tags_per_replicate)) {
    replicate_sample_idx <- experimental_samples[["replicate"]] == rep
    replicate_samples <- experimental_df[, replicate_sample_idx]
    ref_sample_idx <- reference_samples[["replicate"]] == rep
    ref_samples <- as.data.frame(reference_df[, ref_sample_idx])
    rownames(ref_samples) <- rownames(reference_df)
    colnames(ref_samples) <- reference_samplenames[ref_sample_idx]
    observed_tags <- tags_per_replicate[[rep]]
    message("There were ", length(observed_tags), " tags observed in reference replicate ", rep, ".")
    keepers <- rownames(replicate_samples) %in% observed_tags
    replicate_tags <- replicate_samples[keepers, ]
    replicate_tags <- merge(replicate_tags, ref_samples, by = "row.names", all.x = TRUE)
    rownames(replicate_tags) <- replicate_tags[["Row.names"]]
    replicate_tags[["Row.names"]] <- NULL
    na_idx <- is.na(replicate_tags)
    replicate_tags[na_idx] <- 0
    zero_idx <- rowSums(replicate_tags) == 0
    replicate_tags <- replicate_tags[!zero_idx, ]
    message("After filtering for reference tags, there are ", nrow(replicate_tags), " left.")
    replicate_list[[rep]] <- replicate_tags
  }
  return(replicate_list)
}

#' Extract the tags observed in each replicate's reference sample.
#'
#' This ought to be extended to multiple samples/reference.
#'
#' @param meta Metadata matrix.
#' @param reads_per_tag The set of tags from the index list.
gather_tags_per_replicate <- function(meta, reads_per_tag, reference_cutoff = 1) {
  tags_per_replicate <- list()
  reference_sample_idx <- meta[["type"]] == "reference"
  reference_samples <- rownames(meta)[reference_sample_idx]
  replicates_per_reference <- meta[reference_samples, "replicate"]
  names(replicates_per_reference) <- reference_samples
  reference_df <- as.data.frame(reads_per_tag)[, reference_samples]
  rownames(reference_df) <- reads_per_tag[["idx"]]
  replicates <- levels(as.factor(meta[reference_samples, "replicate"]))
  for (rep in 1:length(replicates)) {
    replicate_name <- replicates[rep]
    ref_idx <- meta[["replicate"]] == replicate_name & meta[["type"]] == "reference"
    reference_names <- rownames(meta)[ref_idx]
    if (length(reference_names) == 1) {
      na_idx <- is.na(reference_df[[reference_names]])
      exist <- reference_df[!na_idx, reference_names]
      names(exist) <- rownames(reference_df)[!na_idx]
      sufficient_idx <- exist >= reference_cutoff
      sufficient <- exist[sufficient_idx]
      tags_per_replicate[[replicate_name]] <- names(sufficient)
    } else if (length(reference_names) > 1) {
      sufficient_idx <- rowSums(reference_df[, reference_names], na.rm = TRUE) >= reference_cutoff
      sufficient <- reference_df[sufficient_idx, ]
      tags_per_replicate[[replicate_name]] <- rownames(sufficient)
    } else {
      stop("This replicate does not appear to have any samples associated with it.")
    }
  }
  return(tags_per_replicate)
}
