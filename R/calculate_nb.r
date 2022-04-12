#' Calculate a vector of Nb values from metadata, samples, and references.
#'
#' This function basically retranscribes the original code with new
#' variable names and a few small cleanups.  It takes some metadata,
#' a matrix of sample tags, reference tags, and returns back a vector
#' of the associated Nb values for each sample.
#'
#' There remain quite a few cleanups and improvements which could
#' improve the clarity of this.  Notably, I currently pass the sample
#' names, but those are in the metadata.
#'
#' @param nb_meta Relevant metadata for these samples.
#' @param nb_mtrx Tag matrix for the samples of interest (currently
#'  only calibration samples)
#' @param ref_mtrx Tag matrix for the reference samples.
#' @param nb_samplenames Sample names for the samples of interest
#'  (actually, this provides an easy place to differentiate the
#'  calibration from experimental samples.)
#' @param generations Set the g parameter from the published paper,
#'  this is pretty much always 1.
#' @param zero_filter Remove tags which are not found in the reference
#'  samples?  This should probably be moved out to the various filter
#'  functions.
#' @return Named vector containing the Nb values for the samples of interest.
#' @export
calculate_nb <- function(nb_meta, nb_mtrx, ref_mtrx, nb_samplenames,
                         generations = 1, zero_filter = FALSE) {
  ## Perform a filter for tags not observed in the reference samples
  observed_reference_tags <- rowSums(ref_mtrx, na.rm = TRUE) > 0
  removed_tags <- sum(!observed_reference_tags)
  if (isTRUE(zero_filter)) {
    message("Removed: ", removed_tags, " missing in the reference samples.")
    nb_mtrx <- nb_mtrx[observed_reference_tags, ]
  } else {
    message("Not removing tags which are missing in the reference samples.")
  }

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
  ##  Fhat = 1/k * (Sum from 1 to k) { (fi,s - fi,0)^2 / fi,0 * (1 - fi,0) }

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
  nb_vector <- c()
  for (nb_sample in nb_samplenames) {
    column_vector <- nb_mtrx[, nb_sample]
    column_subset <- column_vector[mean_gt_zero_idx]
    reads_in_sample <- sum(column_subset, na.rm = TRUE)
    sample_frequencies <- column_subset / sum(column_subset, na.rm = TRUE)
    subtracted_frequencies <- initial_frequencies * (1 - initial_frequencies)
    sample_variance <- (sample_frequencies - initial_frequencies) ^ 2
    sample_vs_threshold <- sample_variance / subtracted_frequencies

    Nb_estimate <- generations /
      (mean(sample_vs_threshold, na.rm = TRUE) - 1 / reads_in_sample - 1 / total_reads)
    nb_vector <- c(nb_vector, Nb_estimate)
  }
  names(nb_vector) <- nb_samplenames
  return(nb_vector)
}

#' Given a tag data structure, calculate the Nb values and put them
#' into the metadata.
#'
#' The tag data structure contains the metadata along with the tags.
#' This function therefore passes the relevant pieces to
#' calculate_Nb() and sends the result into the metadata column of
#' choice.  If one specifies the write argument, it will write out a
#' copy of the new metadata.
#'
#' @param tags Data provided by read_csv_tags, read_qiime_otus, or
#'  read_idx_tags.
#' @param generations Set the 'g' parameter from the original paper.
#' @param metadata_column Put the vector of Nb values into this new
#'  column of the metadata.
#' @param zero_filter Filter out tags not observed in the reference
#'  samples?
#' @param write Write the output to this file if not null.
#' @return List containing a backup of the original metadata, new
#'  metadata, the column name with the new Nb values, the vector of Nb
#'  values, and the input tags.
#' @export
nb_from_tags <- function(tags, generations = 1,
                         metadata_column = "Nb", zero_filter=TRUE, write = NULL) {
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

  nb_vector <- calculate_nb(nb_meta, nb_mtrx, ref_mtrx, nb_samplenames,
                            generations=generations, zero_filter=zero_filter)
  for (nb in 1:length(nb_vector)) {
    sample <- names(nb_vector)[nb]
    value <- as.numeric(nb_vector)[nb]
    meta[sample, metadata_column] <- value
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
    "Nb_result" = nb_vector,
    "sample_counts" = tags[["sample_counts"]])
  return(retlist)
}

#' The f(i,s) term in the Nb estimation equation requires frequency estimates.
#'
#' Thus, take the sum of each column and divide every value by it.
#'
#' @param mtrx Matrix to transform.
#' @return Matrix of the same size, containing frequencies.
make_frequency_df <- function(mtrx) {
  freq_mtrx <- mtrx
  for (i in 1:ncol(mtrx)) {
    col_sum <- sum(mtrx[, i], na.rm = TRUE)
    freq_mtrx[, i] <- mtrx[, i] / col_sum
  }
  return(freq_mtrx)
}
