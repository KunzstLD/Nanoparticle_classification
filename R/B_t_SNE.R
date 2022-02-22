# Data prep and t-SNE ----
library(Rtsne)

# Load masses
Masses <- readRDS(file.path(data_cache, "Masses.rds"))

# Scale data
Masses_clust <-
  Masses[!grepl("wafer-|blank-", rownames(Masses)),] %>%
  scale(.) %>%
  as.data.frame(.)

Masses_dist <- Masses_clust %>%
  dist(., method = 'euclidean') %>% 
  as.matrix(.)

# t-SNE
set.seed(42)
tsne_Masses <- Rtsne::Rtsne(X = Masses_dist, perplexity = 20, theta = 0, dims = 2)
tsne_Masses <- as.data.frame(tsne_Masses$Y)

# Create plot ----
# Sample type for plotting
Masses_clust$Sample_type <- paste0("S_", sub("(.+)(\\-)(.+)", "\\1", rownames(Masses_clust)))
Masses_clust$Sample_type <- as.factor(Masses_clust$Sample_type)

sample_color <- c(
  "S_3" = "plum1",
  "S_F3" = "orchid3",
  "S_4" = "cyan",
  "S_F4" = "dodgerblue4",
  "S_9" = "gold",
  "S_F9" = "gold4"
)
ggplot(tsne_Masses, aes(
  x = V1,
  y = V2,
  color = Masses_clust$Sample_type
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
pca_Masses <- prcomp(Masses_clust[-615])
# 54 % variance explained on the first two axis
summary(pca_Masses)
biplot(pca_Masses, scale = TRUE)

library(ggfortify)
autoplot(pca_Masses,
         label = TRUE,
         shape = FALSE)
