#########################################################################
#####  RANDOM FOREST  ######
# Need to run 03_differences_in_intensities.R script prior to this script
#########################################################################

# Add column indicating presence of fulvic acid
Data_RF <-
  cbind(Data, 
        fulvic_acid = as.factor(ifelse(grepl("F", rownames(Data)), 
                                       1, 0)))

# Remove sample 2 (no sorption)
Data_RF <- Data_RF[-which(grepl("F2", rownames(Data))),]

# Remove masses whose intensities decreased
Indexes <- c()
# For selecting the masses increasing for that specific sample
# Increased_Masses <- sub("X","",names(Diffs$S6))
for (i in Increased_Masses) {
  Indexes <- c(Indexes, which(names(Data_RF) == i))
}
Indexes <- c(sort(Indexes), which(names(Data_RF) == "fulvic_acid"))
Data_RF_all <- Data_RF[, Indexes]

# simpler but not ordered
# Dara_RF_all <- Data_RF[, c(Increased_Masses, "fulvic_acid")]

## Random forest analysis function
RF <- function(X,
               seed = 123,
               mtry_step = 1,
               label = NULL) {
  # create training and test data set
  set.seed(seed)
  split <- initial_split(X, prop = .7, strata = "fulvic_acid")
  X_train <- training(split)
  X_test  <- testing(split)
  # number of features
  n_features <- length(setdiff(names(X_train), "fulvic_acid"))
  # Grid for different parameters
  hyper_grid <-
    expand.grid(
      mtry        = seq(2, (length(X) - 1), by = mtry_step),
      node_size   = seq(1, 5, by = 1),
      sample_size = c(0.63, 0.7, 0.8),
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
    hyper_grid$rmse[j] <- sqrt(model$prediction.error)
  }
  
  # top 50 models accoridng to OOB error
  Param_Optimization <-
    hyper_grid[order(hyper_grid$OOB_error), ] %>% head(., 50)
  # use best tuning parameters
  best_set <- hyper_grid[order(hyper_grid$OOB_error),][1,]
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
  names(m_ranger_impurity$variable.importance) <- sub("_", ".",
                                                      sub("Mass_", "",
                                                          names(
                                                            m_ranger_impurity$variable.importance
                                                          )))
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
  names(m_ranger_permutation$variable.importance) <- sub("_", ".",
                                                         sub(
                                                           "Mass_",
                                                           "",
                                                           names(m_ranger_permutation$variable.importance)
                                                         ))
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
    draw_label(
      label = label,
      colour = "darkred",
      size = 23,
      fontface = "bold",
      x = 0.05 ,
      y = 0.96
    )
  
  # plot
  plot_importance <- gridExtra::grid.arrange(p1, p2, nrow = 1)
  # access vip data
  var_imp <-
    vip::vip(m_ranger_impurity,
             num_features = 25,
             geom = NULL)
  RF_masses <- var_imp$data[[1]]
  #### Predict for the test data ####
  pred_class <-
    predict(m_ranger_impurity, X_test[, -which(names(X_test) == "fulvic_acid")])
  # Assess performance on test data
  Confusion_matrix <-
    caret::confusionMatrix(factor(pred_class$predictions),
                           factor(X_test$fulvic_acid))
  
  O <- list(
    Parameters_optimization = Param_Optimization,
    Important_masses = RF_masses,
    Confusion_matrix = Confusion_matrix,
    plot_importance = plot_importance
  )
  return(O)
}

## Random forest analysis for complete dataset
Results_RF_all <- RF(Data_RF_all, mtry_step = 5)

#### Plot Proximity map for complete dataset ####
# Results_RF_all$Parameters_optimization
model1 <-
  randomForest(
    fulvic_acid ~ .,
    data = Data_RF_all,
    ntree = length(Data_RF_all)*10,
    mtry = 27,
    nodesize = 1,
    importance = T,
    proximity = T
  )
MDSplot_modified <-
  function (rf,
            fac,
            k = 2,
            palette = NULL,
            pch = 20,
            ...) {
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
      else
        rainbow(nlevs)
    }
    if (k <= 2) {
      plot(rf.mds$points, col = palette[as.numeric(fac)],
           pch = pch, ...)
      text(
        rf.mds$points,
        labels = rownames(rf.mds$points),
        pos = 3,
        cex = 0.85,
        offset = 0.2
      )
    }
    else {
      pairs(rf.mds$points, col = palette[as.numeric(fac)],
            pch = pch, ...)
    }
    invisible(rf.mds)
  }
MDSplot_modified(model1,
                 Data_RF_all$fulvic_acid,
                 palette = c("blue", "red"),
                 k = 2)

#### Random forest for each individual sunscreen extract ####
Results_RF_Extracts <- list()
Importance_plots <- list()
for (j in 1:length(Diffs_pos)) {
  Indexes <- c()
  #Remove masses whose intensities decreased for one sunscreen extract
  for (i in names(Diffs_pos[[j]])) {
    Indexes <- c(Indexes, which(colnames(Data_RF) == i))
  }
  Indexes <- c(sort(Indexes), which(names(Data_RF) == "fulvic_acid"))
  k <- sub("S", "", names(Diffs_pos)[[j]])
  Data_RF_S <-
    Data_RF[grep(paste(k, "-", sep = ""), rownames(Data_RF)), Indexes]
  Results_RF_Extracts[[j]] <-  RF(Data_RF_S, label = paste0("S", k))
  Importance_plots[[j]] <- Results_RF_Extracts[[j]]$plot_importance
}
names(Results_RF_Extracts) <- names(Diffs_pos)
names(Importance_plots) <- names(Diffs_pos)

# Plot importance plots
ggpubr::ggarrange(plotlist = Importance_plots[1:4],
                  ncol = 2,
                  nrow = 2)
ggpubr::ggarrange(plotlist = Importance_plots[5:8],
                  ncol = 2,
                  nrow = 2)
ggpubr::ggarrange(plotlist = Importance_plots[9],
                  ncol = 2,
                  nrow = 2
)

# Clean up
rm(i, j, Indexes, k)