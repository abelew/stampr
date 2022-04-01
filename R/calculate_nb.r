#' Given the metadata and tag data, calculate Nb values.
#'
#' This is minor cleanup of the original implementation I received from
#' Dr. Lee's collaborators.  I do have the text of the original paper, but it is
#' difficult to interpret; so I fell back on basically retranscribing the
#' original code.
#'
#' @param tags Input tags provided by count_tags.pl and hopefully count_otus.pl
#' @param generations Set the g parameter from the published paper.
#' @export
calculate_nb <- function(tags, generations = 1, metadata_column = "Nb", write = NULL) {
  meta <- tags[["modified_metadata"]]
  meta_backup <- meta
  mtrx <- as.data.frame(tags[["sample_counts"]])
  rownames(mtrx) <- mtrx[["idx"]]
  mtrx[["idx"]] <- NULL
  ## Pull out the non-reference samples for Nb computing
  nb_samples <- meta[["type"]] != "reference"
  nb_samplenames <- rownames(meta)[nb_samples]
  ref_samplenames <- rownames(meta)[!nb_samples]
  ## I am using samplenames as opposed to just the row/column number in case something
  ## got moved around by the user.
  nb_meta <- meta[nb_samplenames, ]
  nb_mtrx <- mtrx[, nb_samplenames]
  ref_mtrx <- mtrx[, ref_samplenames]
  na_idx <- is.na(nb_mtrx)
  nb_mtrx[na_idx] <- 0
  nb_mtrx <- as.matrix(nb_mtrx)
  na_idx <- is.na(ref_mtrx)
  ref_mtrx[na_idx] <- 0
  ref_mtrx <- as.matrix(ref_mtrx)

  ## Perform a filter for tags not observed in the reference samples
  observed_reference_tags <- rowSums(ref_mtrx, na.rm = TRUE) > 0
  removed_tags <- sum(!observed_reference_tags)
  message("Removed: ", removed_tags, " not observed in the reference samples.")
  nb_mtrx <- nb_mtrx[observed_reference_tags, ]

  ## The equation for estimating the bottleneck size is available at:
  ## https://www.nature.com/articles/nmeth.3253 in section 'Bottleneck population size estimate'.
  ## Written out it is:
  ## Nb is similar to Ne which, for each sample, is:
  ##  (The number of generations of competitive growth, g) divided by
  ##    ((Fhat: The mean of the binomial distribution of frequency of reads per tag) - ((1/the sample size) - (1/the inoculum size)))
  ## In turn, Fhat is:
  ##  The sum (for every sample) of: ((each freq./tag - reference freq./tag) squared /
  ##                                  (the reference freq./tag * (1 - the reference freq./tag)))
  ##    divided by the number of samples.
  ##  e.g.:
  ##  Nb ~ Ne = g/(Fhat - 1/S0 - 1/S)
  ##  Fhat = 1/k* (Sum from 1 to k){ (fi,s - fi,0)^2 / fi,0 * (1 - fi,0) }

  ## So, in the equations, fi,0 is the set of mean frequencies which are > 0
  ## from the reference samples.  I will call that term mean_ref_frequencies:
  ref_frequencies <- make_frequency_df(ref_mtrx)
  mean_frequencies <- rowMeans(ref_frequencies, na.rm = TRUE)
  mean_counts <- rowMeans(ref_mtrx, na.rm = TRUE)
  mean_gt_zero_idx <- mean_counts > 0
  mean_frequencies_gt_zero <- mean_frequencies[mean_gt_zero_idx]
  innoculum_frequencies <- rowMeans(ref_frequencies, na.rm = TRUE)
  mean_gt_zero <- mean_counts > 0
  initial_frequencies <- innoculum_frequencies[mean_gt_zero]
  total_reads <- sum(ref_mtrx, na.rm = TRUE)

  ## Set aside a column in the metadata into which to put the Nb estimate
  meta[["Nb"]] <- NA
  for (nb_sample in nb_samplenames) {
    column_vector <- nb_mtrx[, nb_sample]
    column_subset <- column_vector[mean_gt_zero_idx]
    reads_in_sample <- sum(column_subset)
    sample_frequencies <- column_subset / sum(column_subset)
    subtracted_frequencies <- initial_frequencies * (1 - initial_frequencies)
    sample_variance <- (sample_frequencies - initial_frequencies) ^ 2
    sample_vs_threshold <- sample_variance / subtracted_frequencies

    Nb_estimate <- generations /
      (mean(sample_vs_threshold) - 1 / reads_in_sample - 1 / total_reads)
    meta[nb_sample, metadata_column] <- Nb_estimate
  }
  if (!is.null(write)) {
    extension <- tools::file_ext(write)
    if (extension == "xlsx") {
      written <- hpgltools::write_xlsx(data = meta, excel = write)
    } else if (extension == "csv") {
      written <- write.csv(x = meta, file = write)
    } else if (extension == "tsv") {
      writen <- write.tsv(x = meta, file = write)
    } else {
      message("Dunno what to write.")
    }
  }
  retlist <- list(
    "previous_metadata" = meta_backup,
    "modified_metadata" = meta,
    "Nb_column" = metadata_column,
    "sample_counts" = tags[["sample_counts"]])
  return(retlist)
}

#' The f(i,s) term in the Nb estimation equation requires frequency estimates.
#'
#' Thus, take the sum of each column and divide every value by it.
#'
#' @param mtrx Matrix to transform.
make_frequency_df <- function(mtrx) {
  freq_mtrx <- mtrx
  for (i in 1:ncol(mtrx)) {
    col_sum <- sum(mtrx[, i], na.rm = TRUE)
    freq_mtrx[, i] <- mtrx[, i] / col_sum
  }
  return(freq_mtrx)
}
