% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr.R
\name{plot_cor}
\alias{plot_cor}
\title{Plot a heatmap}
\usage{
plot_cor(data, type = c("full", "lower", "upper"), limits_scale = c(-1,
  1), title_legend = "Correlation", palette = c("viridis", "inferno",
  "magma", "plasma"), value = FALSE, color_value = "#ffffff",
  clustering = FALSE, text_size = 12, ...)
}
\arguments{
\item{data}{A matrix or data frame.}

\item{type}{A string to change the type of heatmap. One of "full" (default),
"lower" or "upper".}

\item{limits_scale}{The scale of the fill aesthetic (-1 to 1 by default).}

\item{title_legend}{A string to change the title of the legend.}

\item{palette}{A string indicating the color palette to use. One of
"viridis" (default) "magma", "inferno" or "plasma".}

\item{value}{A boolean (FALSE by default). Should the values should be
visible?}

\item{color_value}{A color string ("#ffffff" by default) for the label
values.}

\item{clustering}{A boolean (FALSE by default). Should the variables be
ordered by a clustering? Use the \code{\link[mclust]{Mclust}} function.}

\item{text_size}{A numeric value to change the size of the labels.}

\item{...}{Parameters to pass to \code{\link[ggplot2]{theme}}.}
}
\value{
A ggplot2 object.
}
\description{
Plot a heatmap of the given data. Usefull to use after
\code{\link{compute_cor}} for instance.
}
\examples{
cor_results <- compute_cor(mtcars)
plot_cor(cor_results)

# change the palette and show the values
plot_cor(cor_results, palette = "inferno", value = TRUE)

# show only the lower triangle and reorder the columns by cluster
plot_cor(cor_results, type = "lower", clustering = TRUE)
}
\seealso{
\code{\link{compute_cor}}
}
