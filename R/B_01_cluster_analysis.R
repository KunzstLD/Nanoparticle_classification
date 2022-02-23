# ___________________________________________________________________________
# DIVISIVE CLUSTERING ----
# ___________________________________________________________________________

# Read in nanoparticle masses
Masses_B_B <- readRDS(file.path(data_cache, "Masses_B.rds"))

# Scale data
Masses_B_clust <-
  Masses_B[!grepl("wafer-|blank-", rownames(Masses_B)),] %>%
  scale(.) %>%
  as.data.frame(.)

# Building of the divisive clusters
Hclust_div <- Masses_B_clust %>%
  dist(., method = 'euclidean') %>%
  cluster::diana(.)

# Plot the dendrogram
factoextra::fviz_dend(Hclust_div, 
                      #k = 12,
                      type = "phylogenic",
                      phylo_layout = "layout.gem", 
                      repel = T,
                      sub = "",
                      cex = 1,
                      lwd = 0.7,
                      ggtheme = theme_void(base_size = 0) 
)

# Different version (hierarchical) with coloring according to the samples
dendro_data <- dendro_data(Hclust_div)
dendro_segments <- dendro_data$segments
setDT(dendro_segments)
dendro_ends <- dendro_segments[yend == 0, ]

# merge label information
dendro_ends[dendro_data$labels,
            `:=`(label = i.label),
            on = "x"]
dendro_ends[, sample_group := sub("(.+)(\\-)(.+)", "\\1", label)]
dendro_ends$sample_group %>% unique()

sample_color <- c(
  "3" = "plum1",
  "F3" = "orchid3",
  "4" = "cyan",
  "F4" = "dodgerblue4",
  "9" = "gold",
  "F9" = "gold4"
)

plot <- ggplot() +
  geom_segment(data = dendro_segments,
               aes(
                 x = x,
                 y = y,
                 xend = xend,
                 yend = yend
               )) +
  geom_segment(data = dendro_ends,
               aes(
                 x = x,
                 y = y,
                 xend = xend,
                 yend = yend,
                 color = sample_group,
                 text = label
               )) +
  scale_color_manual(values = sample_color) +
  scale_y_reverse() +
  coord_flip() +
  labs(x = "", y = "Distance", color = "Sample group") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(family = "Fira", size = 12),
    axis.text.y = element_text(family = "Fira", size = 12),
    legend.title = element_text(family = "Fira", size = 15),
    legend.text = element_text(family = "Fira", size = 12)
  )
ggplotly(plot, tooltip = "text")


# ___________________________________________________________________________
#  K-MEDiOIDS CLUSTERING -----
# ___________________________________________________________________________

# Building K-medioids clusters with different number of clusters 
# and optimizing them using the silhouette width
optsil_obj <- list()
dist_mat <- Masses_B_clust %>% dist(., method = 'euclidean')
for (i in 1:20) {
  K_medioids <- cluster::pam(Masses_B_clust, k = i)
  optsil_obj[[i]] <- optpart::optsil(K_medioids, dist_mat, 100)
  names(optsil_obj)[[i]] <- i
}

# Finding out the optimal number of clusters (k)
Silhouette_widths <-
  lapply(optsil_obj[5:20], function(y)
    summary(cluster::silhouette(y, dist_mat))) %>%
  lapply(., function(y)
    mean(y$clus.avg.widths)) %>%
  unlist(.) %>%
  .[order(.)]
k <- which(Silhouette_widths == max(Silhouette_widths)) %>%
  names(.) %>%
  as.numeric(.)

# Building K-medioids clusters with k (12)
K_medioids <- cluster::pam(Masses_B_clust, k = k)

# Plot clusters
p <- factoextra::fviz_cluster(
  K_medioids,
  data = Masses_B_clust,
  stand = F,
  ellipse = T,
  ellipse.alpha = 0,
  ellipse.type = "convex",
  axes = c(1, 2),
  geom = c("text"),
  repel = T,
  labelsize = 9,
  main = "",
  xlab = F,
  ylab = F,
  ggtheme = theme_void(base_size = 0),
  max.overlaps = 1000
)
p + theme(axis.line = element_blank(),
          legend.position = c(0.2,0.5), legend.direction = "horizontal", legend.title = element_blank(),
          legend.text = element_text(size = 10)
) +
  annotate("text", x = -21, y = -5, 
           label = "k", 
           size = 6)

# Save K-medioids object
saveRDS(object = K_medioids,
        file = file.path(data_cache, "K_mediods_object.rds"))
