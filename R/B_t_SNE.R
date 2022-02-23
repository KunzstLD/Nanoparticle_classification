# Data prep and t-SNE ----
library(Rtsne)

# Load Masses
Masses_B <- readRDS(file.path(data_cache, "Masses_B.rds"))

# Scale data
Masses_B_clust <-
  Masses_B[!grepl("wafer-|blank-", rownames(Masses_B)),] %>%
  scale(.) %>%
  as.data.frame(.)

Masses_B_dist <- Masses_B_clust %>%
  dist(., method = 'euclidean') %>% 
  as.matrix(.)

# t-SNE
set.seed(42)
tsne_Masses_B <- Rtsne::Rtsne(X = Masses_B_dist, perplexity = 20, theta = 0, dims = 2)
tsne_Masses_B <- as.data.frame(tsne_Masses_B$Y)

# Create plot ----
# Sample type for plotting
Masses_B_clust$Sample_type <- paste0("S_", sub("(.+)(\\-)(.+)", "\\1", rownames(Masses_B_clust)))
Masses_B_clust$Sample_type <- as.factor(Masses_B_clust$Sample_type)

sample_color <- c(
  "S_3" = "plum1",
  "S_F3" = "orchid3",
  "S_4" = "cyan",
  "S_F4" = "dodgerblue4",
  "S_9" = "gold",
  "S_F9" = "gold4"
)
ggplot(tsne_Masses_B, aes(
  x = V1,
  y = V2,
  color = Masses_B_clust$Sample_type
)) +
  geom_point(size = 2) +
  scale_color_manual(values = sample_color) +
  labs(color = "Sample type", x = "", y = "") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(family = "Fira", size = 15),
    legend.text = element_text(family = "Fira", size = 12),
  )
ggsave(filename = file.path(data_out, "t_SNE.png"))

# PCA ---
pca_Masses_B <- prcomp(Masses_B_clust[-615])
# 54 % variance explained on the first two axis
summary(pca_Masses_B)
biplot(pca_Masses_B, scale = TRUE)

library(ggfortify)
autoplot(pca_Masses_B,
         label = TRUE,
         shape = FALSE)
