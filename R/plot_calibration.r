#' Plot the calibration curve of the tags observed vs. CFU
#'
#' Currently this is just a simple ggplot scatterplot.  I want to make it more
#' configurable and fun.
#'
#' @param tags Tag data.
#' @param x_column Metadata column for the x axis.
#' @param y_column Metadata column for the y axis.
#' @param color_column Metadata column for colors.
#' @param transform_x Perform a transformation on the x axis?
#' @param transform_y Perform a transformation on the y axis?
#' @export
plot_calibration <- function(tags, x_column = "cfu", y_column = NULL,
                             color_column = "replicate", transform_x = "log10",
                             transform_y = "log10") {
  meta <- tags[["modified_metadata"]]
  calibration_idx <- meta[["type"]] == "calibration"
  cali <- meta[calibration_idx, ]
  if (is.null(y_column)) {
    y_column <- tags[["Nb_column"]]
  }
  if (is.null(meta[[x_column]])) {
    stop("Unable to find the x column: ", x_column, ".")
  }
  if (is.null(meta[[y_column]])) {
    stop("Unable to find the y column: ", y_column, ".")
  }
  x_label <- x_column
  y_label <- y_column

  if (is.null(transform_x)) {
    cali[["x"]] <- cali[[x_column]]
  } else {
    cali[["x"]] <- do.call(what = transform_x, args = list(x = cali[[x_column]]))
    x_label <- glue::glue("{transform_x}({x_column})")
  }
  if (is.null(transform_y)) {
    cali[["y"]] <- cali[[y_column]]
  } else {
    cali[["y"]] <- do.call(what = transform_y, args = list(x = cali[[y_column]]))
    y_label <- glue::glue("{transform_y}({y_column})")
  }
  plot_formula <- y ~ x
  lm_result <- lm(formula = plot_formula, data = cali)
  lm_coef <- coef(lm_result)
  lm_r2 <- summary(lm_result)[["r.squared"]]
  ## glm_result <- glm(formula = plot_formula, data = cali)
  ## glm_coef <- coef(glm_result)
  lm_equation <- glue::glue("y ~ {signif(lm_coef[2], 3)}*x + {signif(lm_coef[1], 3)}, \\
 rsquared: {signif(lm_r2, 3)}")
  a_plot <- ggplot(data = cali,
                   mapping = aes_string(x = "x", y = "y",
                                      fill = color_column, colour = color_column)) +
  scale_shape_manual(values = 21) +
  geom_point(size = 4) +
  ggtitle(lm_equation) + xlab(x_label) + ylab(y_label) +
  geom_rug() +
  geom_smooth(mapping = aes_string(group = "type"), method = "lm", formula = plot_formula)
  ## geom_smooth(mapping = aes_string(group = "type"), method = "glm", formula = plot_formula,
  ##             method.args = list(family = "poisson"))
  retlist <- list(
    "plot" = a_plot,
    "lm" = lm_result)
  return(retlist)
}
