start <- as.POSIXlt(Sys.time())
library(testthat)
library(stampr)
context("001_stampr.R
  1234\n")

## Test that this package works.  I copied a portion of the raw data to inst/extdata.
basedir <- system.file("extdata", package = "stampr")
metadata_file <- file.path(basedir, "all_samples.xlsx")
meta <- hpgltools::extract_metadata(metadata_file)
test_that("We loaded some metadata?", {
  expect_equal(23, nrow(meta))
  expect_equal(9, ncol(meta))
})

## Extract and preprocess the raw reads.
raw_file <- file.path(basedir, "raw_reads.tar")
raw <- utils::untar(tarfile = raw_file)
new_meta <- preprocess_stamp(metadata = meta, raw_column = "rawfastq")
test_that("The metadata has new columns from processing?", {
  expect_equal(11, ncol(new_meta))
})

otu_file <- file.path(basedir, "otus.tar")
raw <- utils::untar(tarfile = otu_file)
otu_mtrx <- read_qiime_otus(new_meta, otu_column = "qiimeotus", trimmed_column = "trimmed_fastq")
test_that("We got a matrix of OTUs from qiime?", {
  expect_equal(24, ncol(otu_mtrx))
  expect_equal(30010, nrow(otu_mtrx))
})



idx_data <- read_tags(new_meta, index_column = "index_table")
idx_meta <- idx_mtrx[["modified_metadata"]]
idx_mtrx <- idx_data[["sample_counts"]]
test_that("We have new data after reading the tags?", {
  expect_equal(23, nrow(idx_meta))
  expect_equal(19, ncol(idx_meta))
  expect_equal(6538, nrow(idx_mtrx))
  expect_equal(24, ncol(idx_mtrx))
})

## Perform an explicit filter using the reference samples and keeping only tags observed in them.
filtered_tags <- filter_given_reference(idx_data)
idx_meta <- filtered_tags[["modified_metadata"]]
idx_mtrx <- filtered_tags[["sample_counts"]]
test_that("We have new data after reference filtering the tags?", {
  expect_equal(23, nrow(idx_meta))
  expect_equal(22, ncol(idx_meta))
  expect_equal(4272, nrow(idx_mtrx))
  expect_equal(24, ncol(idx_mtrx))
})

top500_tags <- filter_topn_tags(filtered_tags)
idx_meta <- top500_tags[["modified_metadata"]]
idx_mtrx <- top500_tags[["sample_counts"]]
test_that("We have new data after top-500 filtering the tags?", {
  expect_equal(23, nrow(idx_meta))
  expect_equal(25, ncol(idx_meta))
  expect_equal(671, nrow(idx_mtrx))
  expect_equal(23, ncol(idx_mtrx))
})

top500_nb <- calculate_nb(top500_tags)
idx_meta <- top500_nb[["modified_metadata"]]
idx_mtrx <- top500_nb[["sample_counts"]]
test_that("We have new data after top-500 filtering the tags?", {
  expect_equal(23, nrow(idx_meta))
  expect_equal(26, ncol(idx_meta))
  expect_equal(671, nrow(idx_mtrx))
  expect_equal(23, ncol(idx_mtrx))
})

top500_plot <- plot_calibration(top500_nb)
test_that("We got a calibration plot?", {
  expect_equal("gg", class(top500_plot)[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 001_stampr.R in ", elapsed,  " seconds."))
