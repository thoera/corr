#' Cramér's V
#'
#' Compute the Cramér's V measure of association.
#'
#' @param x A character vector, matrix or data frame.
#' @param y NULL (default) or a vector if x is a vector.
#' @return A numeric value if \code{x} and \code{y} are vectors or a matrix if
#'   \code{x} is a matrix or a data frame.
#' @seealso \code{\link{compute_cor}}, \code{\link{plot_cor}}
#' @examples
#' compute_cramer_v(tea)
#' compute_cramer_v(tea$price, tea$home)
compute_cramer_v <- function(x, y = NULL) {
  if (is.data.frame(x)) {
    if (sum(vapply(x, function(col) is.factor(col) | is.character(col),
                   logical(1L))) != ncol(x)) {
      stop("'x' must be character or factor")
    }
  }
  if (!(is.data.frame(x) || is.matrix(x)) && is.null(y)) {
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }

  cramer_v <- function(...) {
    test <- chisq.test(...)
    chi2 <- test[["statistic"]]
    n <- sum(test[["observed"]])
    if (test[["method"]] == "Chi-squared test for given probabilities") {
      ind <- which.min(test[["expected"]])
      max_dev <- test[["expected"]]
      max_dev[ind] <- n - max_dev[ind]
      max_chi2 <- sum(max_dev^2L / test[["expected"]])
      v <- sqrt(chi2 / max_chi2)
    }
    else {
      k <- min(dim(test[["observed"]]))
      v <- sqrt(chi2 / (n * (k - 1L)))
    }
    names(v) <- NULL
    return(v)
  }

  if (!(is.data.frame(x) || is.matrix(x))) {
    cramer_results <- suppressWarnings(cramer_v(x, y))
  } else {
    n_cols <- ncol(x)
    names_col <- names(x)
    combinations <- combn(names_col, m = 2L, simplify = FALSE)
    tables <- lapply(combinations, function(col) {
      table(x[, c(col)])
    })
    v <- suppressWarnings(lapply(tables, cramer_v))

    cramer_results <- matrix(1L, nrow = n_cols, ncol = n_cols,
                             dimnames = list(names_col, names_col))
    for (i in seq_along(combinations)) {
      cramer_results[combinations[[i]][1L], combinations[[i]][2L]] <- v[[i]]
      cramer_results[combinations[[i]][2L], combinations[[i]][1L]] <- v[[i]]
    }
  }
  return(cramer_results)
}

#' Correlation and Cramér's V
#'
#' Compute the correlation or the Cramér's V measure of association.
#'
#' @param x A character vector, matrix or data frame.
#' @param y NULL (default) or a vector if x is a vector.
#' @param method a character string indicating which method to use. One of
#'   "pearson" (default), "kendall", "spearman" or "cramer".
#' @param n_jobs numeric, the number of thread(s) to use (1 by default). NB:
#'   n_jobs is not yet implemented.
#' @param ... parameters to pass to \code{\link[stats]{cor}}
#' @return A numeric value if \code{x} and \code{y} are vectors or a matrix if
#'   \code{x} is a matrix or a data frame.
#' @seealso \code{\link{plot_cor}}
#' @examples
#' compute_cor(mtcars)
#' compute_cor(mtcars, method = "spearman")
#' compute_cor(tea$price, tea$home, method = "cramer")
#' @export
compute_cor <- function(x,
                        y = NULL,
                        method = c("pearson", "kendall", "spearman", "cramer"),
                        n_jobs = 1, ...) {
  method <- match.arg(method)

  # if (!is.numeric(n_jobs)) {
  #   stop("'n_jobs' must be numeric", "\n")
  # }

  if (method == "cramer") {
    results <- compute_cramer_v(x = x, y = y)
  } else {
    results <- cor(x = x, y = y, method = method, ...)
  }
  return(results)
}

#' Model-Based Clustering
#'
#' The optimal model according to BIC for EM initialized by hierarchical
#' clustering for parameterized Gaussian mixture models.
#' It uses the \code{\link[mclust]{Mclust}} function.
#' @import mclust
compute_clustering <- function(data) {
  # --
  # Description:
  # Compute the optimal model according to BIC for EM initialized by
  # hierarchical clustering for parameterized Gaussian mixture models.
  # Use the mclust::Mclust() function.
  #
  # Arguments:
  #   data: a numeric vector, matrix or data frame.
  # --

  max_clusters <- ceiling(sqrt(ncol(data)))
  r <- mclust::Mclust(data, G = 2:max_clusters,
                      verbose = FALSE)[["classification"]]
}

#' Plot a heatmap
#'
#' Plot a heatmap of the given data. Usefull to use after
#' \code{\link{compute_cor}} for instance.
#'
#' @param data A matrix or data frame.
#' @param type A string to change the type of heatmap. One of "full" (default),
#'   "lower" or "upper".
#' @param limits_scale The scale of the fill aesthetic (-1 to 1 by default).
#' @param title_legend A string to change the title of the legend.
#' @param palette A string indicating the color palette to use. One of
#'   "viridis" (default) "magma", "inferno" or "plasma".
#' @param value A boolean (FALSE by default). Should the values should be
#'   visible?
#' @param color_value A color string ("#ffffff" by default) for the label
#'   values.
#' @param clustering A boolean (FALSE by default). Should the variables be
#'   ordered by a clustering?
#' @param text_size A numeric value to change the size of the labels.
#' @param ... Parameters to pass to \code{\link[ggplot2]{theme}}.
#' @return A ggplot2 object.
#' @seealso \code{\link{compute_cor}}
#' @examples
#' cor_results <- compute_cor(mtcars)
#' plot_cor(cor_results)
#'
#' # change the palette and show the values
#' plot_cor(cor_results, palette = "inferno", value = TRUE)
#'
#' # show only the lower triangle and reorder the columns by cluster
#' plot_cor(cor_results, type = "lower", clustering = TRUE)
#' @export
plot_cor <- function(data,
                     type = c("full", "lower", "upper"),
                     limits_scale = c(-1, 1),
                     title_legend = "Correlation",
                     palette = c("viridis", "inferno", "magma", "plasma"),
                     value = FALSE,
                     color_value = "#ffffff",
                     clustering = FALSE,
                     text_size = 14, ...) {
  if (!(is.matrix(data) | is.data.frame(data))) {
    stop("'data' must be a matrix or a data frame")
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  type <- match.arg(type)
  palette <- match.arg(palette)
  if (!is.logical(value)) {
    stop("'value' must be logical")
  }

  # if clustering == TRUE, compute the clustering to find how to order
  # the variables
  if (clustering == TRUE) {
    clusters <- compute_clustering(data)
    clustering_order <- names(sort(clusters))

    # reorder the rows and columns
    data <- data[clustering_order, clustering_order]
  }

  # if type == "lower" or type == "upper" remove one side of the matrix
  if (type == "lower") {
    data[upper.tri(data)] <- NA
  }
  if (type == "upper") {
    data[lower.tri(data)] <- NA
  }

  # reshape the data to long format
  cor_ggplot <- cbind(expand.grid(dimnames(data)), value = as.vector(data))
  cor_ggplot <- cor_ggplot[!is.na(cor_ggplot[["value"]]), ]

  # a ggplot2 theme
  theme_heatmap <- function(text_size = text_size) {
    ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(size = text_size),
                     axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1,
                                                         vjust = 0.5),
                     legend.key.height = grid::unit(0.075, "npc"),
                     legend.key.width = grid::unit(0.035, "npc"),
                     legend.text = ggplot2::element_text(size = 12),
                     legend.title = ggplot2::element_text(size = 14))
  }

  g <- ggplot2::ggplot(data = cor_ggplot,
                       ggplot2::aes(x = Var1, y = Var2, fill = value,
                                    text = paste0(
                                      "x: ", Var1, "\n",
                                      "y: ", Var2, "\n",
                                      "Value: ", round(value, 2L))
                       )) +
    ggplot2::geom_tile(color = "#ffffff") +
    viridis::scale_fill_viridis(limits = limits_scale,
                                option = palette,
                                name = title_legend) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "", x = "", y = "") +
    theme_heatmap(text_size = text_size) +
    ggplot2::theme(...)
  if (type %in% c("full", "lower")) {
    g <- g +
      ggplot2::scale_x_discrete(limits = rev(rownames(data)))
  }
  if (type == "upper") {
    g <- g +
      ggplot2::scale_x_discrete(position = "top",
                                limits = rev(rownames(data))) +
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::theme(legend.position = "left",
                     axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 0,
                                                         vjust = 0.5))
  }
  if (value == TRUE) {
    g <- g + ggplot2::geom_text(ggplot2::aes(label = round(value, 2L)),
                                color = color_value, size = text_size / 3L)
  }
  return(g)
}
