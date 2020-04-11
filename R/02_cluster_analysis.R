###########################
### DIVISIVE CLUSTERING ###
###########################

# Scale data
Data_Clust <- scale(Data)
Data_Clust <- as.data.frame(Data_Clust)

# Building of the divisive clusters
dist_mat <- dist(Data_Clust, method = 'euclidean')
Hclust_div <- cluster::diana(dist_mat)

# Plot the dendrogram
factoextra::fviz_dend(Hclust_div, 
                      k = 24,
                      type = "phylogenic",
                      phylo_layout = "layout.gem", 
                      repel = T,
                      sub = "",
                      cex = 1.3,
                      lwd = 0.7,
                      ggtheme = theme_void(base_size = 0) 
)

# Clean up
rm(Hclust_div)

# TODO: Outliers in Clustering?

###########################
#  K-MEDiOIDS CLUSTERING  #
###########################

# Building K-medioids clusters with different number of clusters 
# and optimizing them using silhouette wigth
optsil_obj <- list()
for(i in 1:20) {
  K_medioids <- pam(Data_Clust, k = i)
  optsil_obj[[i]] <- optpart::optsil(K_medioids, dist_mat, 100)
  names(optsil_obj)[[i]] <- i
}

# Finding out the optimal number of clusters (k)
Silhouette_widths <-
  lapply(optsil_obj[5:20], function(y)
    summary(silhouette(y, dist_mat))) %>%
  lapply(., function(y)
    mean(y$clus.avg.widths)) %>%
  unlist(.) %>%
  .[order(.)]
k <- which(Silhouette_widths == max(Silhouette_widths)) %>%
  names(.) %>%
  as.numeric(.)

# Building K-medioids clusters with k
K_medioids <- pam(Data_Clust, k=k)

# Plot clusters
p <- fviz_cluster(K_medioids, data = Data_Clust, 
                  stand = F, 
                  ellipse = T, ellipse.alpha = 0, ellipse.type = "convex",
                  axes = c(1,2),
                  geom = c("text"),
                  repel = T,
                  labelsize = 12,
                  main = "",
                  xlab = F, ylab = F,
                  ggtheme = theme_void(base_size = 0)
)
p + theme(axis.line = element_blank(),
          legend.position = c(0.2,0.45), legend.direction = "horizontal", legend.title = element_blank(),
          legend.text = element_text(size = 10)
) +
  annotate("text", x = -18.4, y = -5.8, 
           label = "k", 
           size = 6)

# Clean up
rm(i, dist_mat,Silhouette_widths, p,k,optsil_obj,K_medioids)