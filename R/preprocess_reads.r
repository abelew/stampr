#' Preprocess raw reads using cutadapt and my little perl script.
#'
#' @param metadata Sample metadata.
#' @param column Metadata column containing the locations of the raw reads.
#' @export
preprocess_stamp <- function(metadata, raw_column = "raw_fastq",
                            trimmed_column = "trimmed_fastq",
                            index_column = "index_table", output = NULL) {
  if ("character" %in% class(metadata)) {
    metadata <- hpgltools::extract_metadata(metadata)
  }

  do_trim <- TRUE
  do_count <- TRUE
  trimmed_files <- NULL
  tag_files <- NULL
  raw_files <- metadata[[raw_column]]
  if (is.null(raw_files)) {
    stop("This requires a column of the input raw data.")
  }
  raw_exist <- file.exists(raw_files)
  if (sum(raw_exist) < length(raw_files)) {
    stop("Some of the raw files appear to be missing.")
  }

  ## Start out checking if this is already finished.
  if (!is.null(metadata[[trimmed_column]])) {
    trimmed_files <- metadata[[trimmed_column]]
    trimmed_exist <- file.exists(trimmed_files)
    if (sum(trimmed_exist) == nrow(metadata)) {
      do_trim <- FALSE
    }
  }

  if (!is.null(metadata[[index_column]])) {
    index_files <- metadata[[index_column]]
    index_exist <- file.exists(index_files)
    if (sum(index_exist) == nrow(metadata)) {
      do_count <- FALSE
    }
  }

  if (isTRUE(do_trim)) {
    for (f in 1:length(raw_files)) {
      file <- raw_files[f]
      trimmed <- metadata[f, trimmed_column]
      dname <- dirname(file)
      bname <- basename(file)
      tbase <- gsub(x = bname, pattern = "\\.fastq.*$", replacement = "")
      trimmed <- glue::glue("{tbase}_trimmed.fastq")
      cmd <- glue::glue("less {file} | cutadapt - -a GTCATAGCTGTTT -e 0.1 -n 3 -m 8 -M 42 --too-short-output={dname}/tooshort.fastq.gz --too-long-output={dname}/toolong.fastq.gz --untrimmed-output={dname}/untrimmed.fastq.gz -o {dname}/{trimmed} 2>&1 1>{dname}/{tbase}_cutadapt.log")
      message("Running: ", cmd, ".")
      trim_result <- system2("/bin/bash", args = c("-c", shQuote(cmd)), stdout = TRUE)
      cmd <- glue::glue("xz -9e -f {dname}/{trimmed}")
      xz_result <- system2("/bin/bash", args = c("-c", shQuote(cmd)), stdout = TRUE)
      metadata[f, trimmed_column] <- glue::glue("{dname}/{trimmed}.xz")
    }  ## End iterating over every raw file.
  }  ## End trimming

  trimmed_files <- metadata[[trimmed_column]]
  if (isTRUE(do_count)) {
    for (f in 1:length(trimmed_files)) {
      file <- trimmed_files[f]
      dname <- dirname(file)
      bname <- basename(file)
      count <- gsub(x = bname, pattern = "\\.fastq.*$", replacement = "")
      counted <- glue::glue("{dname}/{count}_tags.txt")
      counter <- system.file("extdata/count_tags.pl", package = "stampr")

      cmd <- glue::glue("{counter} {file} {counted}")
      message("Running: ", cmd, ".")
      count_result <- system2("/bin/bash", args = c("-c", shQuote(cmd)), stdout = TRUE)
      cmd <- glue::glue("xz -9e -f {counted}")
      xz_result <- system2("/bin/bash", args = c("-c", shQuote(cmd)), stdout = TRUE)
      counted <- paste0(counted, ".xz")
      metadata[f, index_column] <- counted
    }  ## End iterating over every trimmed file.
  }  ## End counting
  return(metadata)
}
