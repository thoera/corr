# corr

## Overview

corr aims to provide a consistent interface to compute several types of correlation coefficient and Cramér's V measure of association. It offers an easy way to visualize the results with heatmaps.

corr fits naturally in the `tidyverse` and is "pipe-friendly".

## Installation

To install the package, simply run the following from an R console:

```r
# install.packages("devtools")
devtools::install_github("thoera/corr")
```

## Usage

```r
library("corr")

# compute the correlation
results <- compute_cor(mtcars, method = "pearson")

# view the results 
plot_cor(results, value = TRUE)

# corr is "pipe-friendly"
library("tidyverse")

starwars %>%
  select_if(is.character) %>%
  select(-name) %>%
  filter(complete.cases(.)) %>%
  compute_cor(method = "cramer") %>%
  plot_cor(type = "lower",
           value = TRUE,
           limits_scale = c(0, 1),
           title_legend = "Cramér's V:")
```

## Example

With the `mtcars` dataset:

![examples.png](/README_examples/examples.png?raw=true)
