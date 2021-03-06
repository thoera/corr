cor_results <- compute_cor(mtcars, method = "pearson")

png("README_examples/examples.png", width = 1920, height = 1080, res = 150)
plot_cor(cor_results,
         clustering = TRUE,
         title_legend = "Correlation (Pearson):")
dev.off()
