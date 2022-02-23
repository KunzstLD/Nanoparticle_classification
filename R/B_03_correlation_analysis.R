# ___________________________________________________________________________
#  CORRELATION ANALYSIS  ----
# ___________________________________________________________________________

# Read in nanoparticle masses and Diffs file (significant differences between certain masses)
Masses_B <- readRDS(file.path(data_cache, "Masses_B.rds"))
Diffs <- readRDS(file.path(data_cache, "Diffs.rds"))

Output_cormat <- list()
for (i in names(Diffs)[c(2, 3, 4)]) {
  # Determine the correlation matrix for the increased Masses_B of the sample in the whole data set
  Cor_matrix <-
    Masses_B[, Diffs[[i]]$`Increased Masses_B`] %>% cor(., method = "pearson") %>%
    as.data.frame(.)
  colnames(Cor_matrix) <-
    Cor_matrix %>% colnames(.) %>% sub("Mass_", "", .) %>% sub("_", "", .)
  rownames(Cor_matrix) <-
    Cor_matrix %>% rownames(.) %>% sub("Mass_", "", .) %>% sub("_", "", .)
  Output_cormat[[i]] <- Cor_matrix
}
Output_cormat[[1]]

# Correlation plot for specific Masses_B
# Choose only Masses_B that have high positive correlations (> 0.9), 
# here at least 50 out of 281 (for F3-3)
high <- apply(Output_cormat$`F3-3`, 2, function(y)
  sum(y > 0.9)) %>%
  .[. >= 50] %>%
  names(.)

f3_cor_matrix <- as.matrix(Output_cormat$`F3-3`)

png(
  filename = file.path(data_out, "example_corr_plot.png"),
)
corrplot::corrplot(
  as.matrix(f3_cor_matrix[high, high]),
  method = "color",
  type = "lower",
  tl.col = "black",
  tl.cex = .7,
  order = "hclust"
)
dev.off()