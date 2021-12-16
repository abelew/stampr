#' Create data structures similar to read_idx using qiime otus.
#'
#' There are some problems with this implementation still, primarily because the
#' version of qiime does not helpfully give the tag sequences.  I have found
#' some ways around this, but since I don't really like that implementation I
#' haven't finished it yet.
#'
#' @param metadata Sample metadata.
#' @param column Metadata column containing the qiime output.
#' @param output Write the matrix to this file, if provided.
#' @export
read_qiime_otus <- function(metadata, otu_column = "qiime_otus", xref_sequence = FALSE,
                            trimmed_column = "trimmed_reads", output = NULL) {
  modified <- metadata
  modified[["qiime_reads"]] <- 0
  modified[["qiime_observed"]] <- 0
  modified[["qiime_max_reads_per_tag"]] <- 0
  otu_files <- metadata[[otu_column]]
  if (is.null(otu_files)) {
    stop("The input file column does not include the qiime outputs.")
  }
  trimmed_files <- metadata[[trimmed_column]]
  if (is.null(trimmed_files)) {
    stop("The input file column does not include the trimmed reads.")
  }
  count_names <- gsub(x = basename(otu_files), pattern = "_otus.*", replacement = "")
  seq_names <- paste0(count_names, "_seq")
  final_dt <- data.table::data.table()
  for (f in 1:length(otu_files)) {
    otu_file <- otu_files[f]
    trimmed_file <- trimmed_files[f]
    message("Starting to read: ", otu_file, ".")
    count_name <- gsub(x = basename(otu_file), pattern = "_otus.*", replacement = "")
    sequence_name <- paste0(count_name, "_seq")
    column_names <- c(sequence_name, count_name)
    interim_column_names <- c("otu_name", "number_reads", "represent_read")
    ## I am using readLines() instead of something like read.table() because the qiime output
    ## is not in a standard tsv/etc format, but instead starts with 1 column
    ## containing the otu name followed by a variable number of columns, 1 for
    ## every read comprising the otu.  The purpose of the tally script is to count
    ## these and provide a 2 column matrix / sample containing the name and number
    ## of reads.
    lines <- readLines(otu_file)
    pair <- data.table::as.data.table(matrix(nrow = length(lines), ncol = 3))
    pair[[1]] <- as.character(pair[[1]])
    pair[[2]] <- as.numeric(pair[[2]])
    pair[[3]] <- as.character(pair[[3]])
    colnames(pair) <- interim_column_names
    for (l in 1:length(lines)) {
      line <- lines[l]
      split <- strsplit(x = line, split = "\t")[[1]]
      pair[l, 1] <- split[1]
      pair[l, 3] <- split[2]
      ## The -1 is because the first element is the name of the OTU, thus the
      ## length-1 is the number of reads comprising the OTU.
      pair[l, 2] <- length(split) - 1
    }
    modified[f, "qiime_reads"] <- sum(pair[, 2])
    modified[f, "qiime_observed"] <- nrow(pair)
    modified[f, "qiime_max_reads_per_tag"] <- max(pair[, 2])

    if (isTRUE(xref_sequence)) {
      message("Cross referencing otus with tag IDs.")
      cmd <- as.character(glue::glue("seqtk subseq -t <(less {trimmed_file}) <(awk '{{print $2}}' <(less {otu_file}))"))
      sequence_ids <- system2("/bin/bash", args = c("-c", shQuote(cmd)), stdout = TRUE)
      sequence_df <- readr::read_tsv(file = sequence_ids, col_names = c("ID", "number", "sequence"))
      pair <- merge(pair, sequence_df, by.x = "represent_read", by.y = "ID")
      pair <- pair[, c("sequence", "number_reads")]
      colnames(pair) <- c("sequence", count_name)
    } else {
      pair <- pair[, c("otu_name", "number_reads")]
      colnames(pair) <- c("sequence", count_name)
    }

    ## Make a big matrix from all samples
    if (f == 1) {
      final_dt <- pair
    } else {
      final_dt <- merge(final_dt, pair, by = "sequence", all = TRUE, allow.cartesian = TRUE)
      ##print(dim(final_dt))
    }
  }  ## Finished iterating over every otu/trimmed file.
  if (!is.null(output)) {
    write.csv(x = final_dt, file = output)
  }
  retlist <- list(
    "modified_metadata" = modified,
    "sample_counts" = final_dt)
  return(retlist)
}
