#' Make a boxplot of read density
#'
#' @param retlist from read_idx or read_qiime
#' @export
plot_read_density <- function(retlist) {
  long_samples <- data.frame()
  if (is.null(retlist[["long_samples"]])) {
    ## Then this came from qiime_otus.
    long_samples <- reshape2::melt(retlist[["sample_counts"]])
    colnames(long_samples) <- c("sequence", "sampleid", "number")
    met <- retlist[["modified_metadata"]][, c("sampleid", "replicate")]
    long_samples <- merge(long_samples, met, by = "sampleid")
  } else {
    long_samples <- retlist[["long_samples"]]
  }
  plt <- ggplot(data = long_samples,
                mapping = aes_string(fill = "sampleid", x = "sampleid", y = "number")) +
  geom_boxplot() +
  scale_y_log10()
  return(plt)
}
