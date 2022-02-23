# ___________________________________________________________________________
# Analysis A 
# Presented in one script, might be split up in several scripts later
# ___________________________________________________________________________
Masses_A <- readRDS(file.path(data_cache, "Masses_A.rds"))

# ___________________________________________________________________________
# DIVISIVE CLUSTERING ----
# ___________________________________________________________________________

# Scale Masses_A
Masses_A_Clust <- scale(Masses_A)
Masses_A_Clust <- as.data.frame(Masses_A_Clust)

# Building of the divisive clusters
dist_mat <- dist(Masses_A_Clust, method = 'euclidean')
Hclust_div <- cluster::diana(dist_mat)

# Plot the dendrogram
factoextra::fviz_dend(Hclust_div, 
                      #k = 24,
                      type = "phylogenic",
                      phylo_layout = "layout.gem", 
                      repel = T,
                      sub = "",
                      cex = 1.3,
                      lwd = 0.7,
                      ggtheme = theme_void(base_size = 0) 
)

# ___________________________________________________________________________
#  K-MEDiOIDS CLUSTERING  ----
# ___________________________________________________________________________

# Building K-medioids clusters with different number of clusters and optimizing them using silhouette wigth
optsil_obj <- list()
for(i in 1:20) {
  K_medioids <- cluster::pam(Masses_A_Clust, k = i)
  optsil_obj[[i]] <- optpart::optsil(K_medioids, dist_mat, 100)
  names(optsil_obj)[[i]] <- i
}

# Finding out the optimal number of clusters (k)
Silhouette_widths <- lapply(optsil_obj[5:20], function(y) summary(silhouette(y, dist_mat))) %>%
  lapply(., function(y) mean(y$clus.avg.widths)) %>%
  unlist(.) %>%
  .[order(.)]
k <- which(Silhouette_widths == max(Silhouette_widths)) %>%
  names(.) %>%
  as.numeric(.)

# Building K-medioids clusters with k
K_medioids <- pam(Masses_A_Clust, k=k)

# Plot clusters
p <- factoextra::fviz_cluster(K_medioids, data = Masses_A_Clust, 
                  stand = F, 
                  ellipse = T, ellipse.alpha = 0, ellipse.type = "convex",
                  axes = c(1,2),
                  geom = c("text"),
                  repel = T,
                  labelsize = 9,
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


# ___________________________________________________________________________
#  DIFFERENCES IN INTENSITIES  ----
# ___________________________________________________________________________

# Determine and list significant differences between exposed samples (marked with "F") and non-exposed samples
Diffs <- list()
Increased_Masses <- c()

for (i in c(1:10)){
  # Select the Masses_A corresponding to fulvic acid exposure and non-exposed
  Ful   <- Masses_A[grep(rownames(Masses_A), pattern = paste0("F",i,"-[1-9]")),]
  NoFul <- Masses_A[grep(rownames(Masses_A), pattern = paste0("^",i,"-[0-9]")),]
  # Calculate the differences of the medians and the mean values of the non-exposed samples
  Diff <-  apply(Ful, FUN = median, MARGIN = 2) - apply(NoFul, FUN = median, MARGIN = 2)
  Means <- apply(NoFul, FUN = mean , MARGIN = 2)
  # Select the significant (>95%) differences and normalize them
  for (j in 1:length(Diff)){
    if (t.test(Ful[[j]], NoFul[[j]])$p.value > 0.05){Diff[[j]] <- NA}
    else {Diff[[j]] <- Diff[[j]]/Means[[j]]
    # Make a list of masses whose intensites incresed significantly 
    if (Diff[[j]] > 0){Increased_Masses <- c(Increased_Masses, names(Diff)[[j]])}
    }
  }
  Diffs[[i]] <- Diff
}
names(Diffs) <- paste0("S",1:10)
# Remove duplicates in the list of masses whose intensites increased
Increased_Masses <- Increased_Masses[!duplicated(Increased_Masses)]

# Preparing for ploting
Diff_plot <- data.frame()
for (i in c(1:length(Diffs))) {
  A <- data.frame(
    diff = Diffs[[i]],
    mass = as.numeric(sub("_", ".", sub("Mass_", "", names(Diffs[[i]])))),
    sample = as.factor(rep(names(Diffs)[[i]], length(Diffs[[i]])))
  )
  Diff_plot <- rbind(Diff_plot, A)
}

# Plot
ggplot(Diff_plot, aes(x=mass, y=diff))+
  geom_point(size = 0.8, shape = 21)+
  facet_wrap(~sample, 
             scale = "free_y", 
             ncol = 3) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11,face="bold"),
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"))+
  labs(title = element_blank(),
       x="m/z",
       y="Normalized intensity differences", 
       linetype=NULL)

# Create a list of increased intensities for each separate sample
Diffs_pos <- Diffs
for (i in 1:length(Diffs)){Diffs_pos[[i]] <- Diffs[[i]][which(Diffs[[i]] > 0)]}

# Determine the common increasing masses for each sunscreen extract:
Common_masses <- as.data.frame(matrix(0, ncol = length(Diffs_pos), nrow = length(Diffs_pos)))
for (i in 1:length(Diffs_pos)){
  for (j in 1:length(Diffs_pos)){
    Common_masses[j,i] <- length(intersect(names(Diffs_pos[[i]]), names(Diffs_pos[[j]])))
  }
}
# Remove values above the diagonal and give names to the columns/rows
Common_masses[upper.tri(Common_masses)] <- ""
colnames(Common_masses) <- names(Diffs_pos)
rownames(Common_masses) <- names(Diffs_pos)

# Export
write.csv(file = file.path(Masses_A_out, "Common masses.csv"), x = Common_masses)


# ___________________________________________________________________________
#  CORRELATION ANALYSIS  ----
# ___________________________________________________________________________

# Remove S2 (no increased masses)
Diffs_pos$S2 <- NULL

# Build a list of network plots
list_plot <- list()
for (i in c(1:length(Diffs_pos))){
  # Determine the correlation matrix for the increased masses of the sample in the whole Masses_A set
  Cor.matrix <- cor(Masses_A[,names(Diffs_pos[[i]])], method = "pearson")
  colnames(Cor.matrix) <- sub("_",".",sub("Mass_","",colnames(Cor.matrix)))
  rownames(Cor.matrix) <- sub("_",".",sub("Mass_","",rownames(Cor.matrix)))
  # Save the network plot
  list_plot[[i]] <- network_plot(Cor.matrix,min_cor = 0.7,repel = T)+
    theme(legend.text = element_text(size = 16, face = "bold"),
          legend.key.size = unit(1.3, "cm"))
}
ggpubr::ggarrange(plotlist = list_plot,
                  ncol = 2, nrow = 5,
                  font.label = list(size = 26, face = "bold", color = "darkred"),
                  labels = names(Diffs_pos),
                  legend = "top",
                  common.legend = T)


# ___________________________________________________________________________
#  RANDOM FOREST  ----
# ___________________________________________________________________________

# Add column indicating presence of fulvic acid
Masses_A_RF <- cbind(Masses_A, fulvic_acid = as.factor(ifelse(grepl("F", rownames(Masses_A)), 1, 0)))
# Remove sample 2 (no sorption)
Masses_A_RF <- Masses_A_RF[-which(grepl("F2", rownames(Masses_A))),]

# Remove masses whose intensities decreased
Indexes <- c()
# For selecting the masses increasing for that specific sample
# Increased_Masses <- sub("X","",names(Diffs$S6))
for (i in Increased_Masses) {
  Indexes <- c(Indexes, which(colnames(Masses_A_RF) == i))
}
Indexes <- c(sort(Indexes), which(names(Masses_A_RF) == "fulvic_acid"))
Masses_A_RF_all <- Masses_A_RF[, Indexes]

## Random forest analysis function
RF <- function(X, seed = 123, mtry_step = 1, label = NULL){
  # create training and test Masses_A set
  set.seed(seed)
  split <- initial_split(X, prop = .7, strata = "fulvic_acid")
  X_train <- training(split)
  X_test  <- testing(split)
  # number of features
  n_features <- length(setdiff(names(X_train), "fulvic_acid"))
  # Grid for different parameters
  hyper_grid <- expand.grid(mtry        = seq(1, (length(X)-1), by = mtry_step),
                            node_size   = seq(1, 5, by = 1),
                            sample_size = c(0.63,0.7,0.8),
                            rmse = NA,
                            OOB_error = NA
  )
  # Excute RF with grid search
  for (j in 1:nrow(hyper_grid)) {
    # train model
    model <- ranger(
      formula         = fulvic_acid ~ .,
      data            = X_train,
      seed            = 123,
      verbose         = FALSE,
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$node_size[j],
      sample.fraction = hyper_grid$sample_size[j]
    )
    # add OOB error to grid
    hyper_grid$OOB_error[j] <- model$prediction.error
  }
  
  # top 50 models accoridng to OOB error
  Param_Optimization <- hyper_grid[order(hyper_grid$OOB_error),] %>% head(., 50)
  # use best tuning parameters
  best_set <- hyper_grid[order(hyper_grid$OOB_error), ][1, ]
  # re-run model with impurity-based variable importance
  m_ranger_impurity <- ranger(
    formula         = fulvic_acid ~ .,
    data            = X_train,
    num.trees       = n_features * 10,
    mtry            = best_set$mtry,
    min.node.size   = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance      = 'impurity',
    verbose         = FALSE,
    seed            = 123
  )
  
  # re-run model with permutation-based variable importance
  m_ranger_permutation <- ranger(
    formula         = fulvic_acid ~ .,
    data            = X_train,
    num.trees       = n_features * 10,
    mtry            = best_set$mtry,
    min.node.size   = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance      = 'permutation',
    verbose         = FALSE,
    seed            = 123
  )
  # most important variables according to impurtiy-based VI in distinguishing TPGs
  names(m_ranger_impurity$variable.importance) <- sub("_",".",
                                                      sub("Mass_","",
                                                          names(m_ranger_impurity$variable.importance)))
  p1 <- vip::vip(m_ranger_impurity,
                 num_features = 25,
                 geom = "point") +
    ggtitle(paste("Impurity-based")) +
    geom_point(size = 3) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 15),
      axis.text.x = element_text(family = "Fira", size = 11),
      axis.text.y = element_text(family = "Fira", size = 11)
    ) +
    xlab("Mass in m/z")
  
  # most important variables according to permutation-based VI in distinguishing TPGs
  names(m_ranger_permutation$variable.importance) <- sub("_",".",
                                                         sub("Mass_","",
                                                             names(m_ranger_permutation$variable.importance)))
  p2 <- vip::vip(m_ranger_permutation,
                 num_features = 25,
                 geom = "point") +
    ggtitle(paste("Permutation-based")) +
    geom_point(size = 3) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 15),
      axis.text.x = element_text(family = "Fira", size = 11),
      axis.text.y = element_text(family = "Fira", size = 11)
    )
  p2 <- ggdraw(p2) + 
    draw_label(label = label, colour = "darkred", size = 23, fontface = "bold", x = 0.05 , y = 0.96)
  
  # plot
  plot_importance <- gridExtra::grid.arrange(p1, p2, nrow = 1)
  # access vip Masses_A 
  var_imp <- vip::vip(m_ranger_impurity, num_features = 25, geom = NULL)
  RF_masses <- var_imp$data[[1]]
  #### Predict for the test Masses_A ####
  pred_class <- predict(m_ranger_impurity, X_test[, -which(names(X_test) == "fulvic_acid")])
  # Assess performance on test Masses_A
  Confusion_matrix <- caret::confusionMatrix(factor(pred_class$predictions), 
                                             factor(X_test$fulvic_acid))
  
  O <- list(Parameters_optimization = Param_Optimization, 
            Important_masses = RF_masses, 
            Confusion_matrix = Confusion_matrix,
            plot_importance = plot_importance) 
  return(O)
}

## Random forest analysis for complete Masses_Aset
Results_RF_all <- RF(Masses_A_RF_all, mtry_step = 5)

# Plot Proximity map for complete Masses_Aset
library(randomForest)
model1 <- randomForest(fulvic_acid ~ ., Data = Masses_A_RF_all, ntree = 1200, mtry = 21, nodesize = 1, importance = T, proximity = T)
MDSplot_modified <- function (rf, fac, k = 2, palette = NULL, pch = 20, ...) {
  if (!inherits(rf, "randomForest")) 
    stop(deparse(substitute(rf)), " must be a randomForest object")
  if (is.null(rf$proximity)) 
    stop(deparse(substitute(rf)), " does not contain a proximity matrix")
  op <- par(pty = "s")
  on.exit(par(op))
  rf.mds <- stats::cmdscale(1 - rf$proximity, eig = TRUE, 
                            k = k)
  colnames(rf.mds$points) <- paste0("Dimension ", 1:k)
  nlevs <- nlevels(fac)
  if (is.null(palette)) {
    palette <- if (requireNamespace("RColorBrewer", quietly = TRUE) && 
                   nlevs < 12) 
      RColorBrewer::brewer.pal(nlevs, "Set1")
    else rainbow(nlevs)
  }
  if (k <= 2) {
    plot(rf.mds$points, col = palette[as.numeric(fac)], 
         pch = pch, ...)
    text(rf.mds$points, labels = rownames(rf.mds$points), pos = 3, cex = 0.85, offset = 0.2)
  }
  else {
    pairs(rf.mds$points, col = palette[as.numeric(fac)], 
          pch = pch, ...)
  }
  invisible(rf.mds)
}
MDSplot_modified(model1, Masses_A_RF_all$fulvic_acid, palette=c("blue","red"), k = 2)

## Random forest for each individual sunscreen extract
Results_RF_Extracts <- list()
Importance_plots <- list()
for (j in 1:length(Diffs_pos)) {
  Indexes <- c()
  # Remove masses whose intensities decreased for one sunscreen extract
  for (i in names(Diffs_pos[[j]])) {
    Indexes <- c(Indexes, which(colnames(Masses_A_RF) == i))
  }
  Indexes <- c(sort(Indexes), which(names(Masses_A_RF) == "fulvic_acid"))
  k <- sub("S", "", names(Diffs_pos)[[j]])
  Masses_A_RF_S <- Masses_A_RF[grep(paste(k, "-", sep = ""), rownames(Masses_A_RF)), Indexes]
  Results_RF_Extracts[[j]] <- RF(Masses_A_RF_S, label = paste0("S", k))
  Importance_plots[[j]] <- Results_RF_Extracts[[j]]$plot_importance
}
names(Results_RF_Extracts) <- names(Diffs_pos)
names(Importance_plots) <- names(Diffs_pos)

# Plot importance plots
ggpubr::ggarrange(plotlist = Importance_plots[1:4],
                  ncol = 2, nrow = 2
)
ggpubr::ggarrange(plotlist = Importance_plots[5:8],
                  ncol = 2, nrow = 2
)
ggpubr::ggarrange(plotlist = Importance_plots[9],
                  ncol = 2, nrow = 2
)