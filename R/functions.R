# ___________________________________________________________________________
# Functions used for data processing & analysis of nanoparticle masses ----
# ___________________________________________________________________________

## Data processing ----
# Function to import and prepare the data to be combined including normalization
Import <- function(X, Normalization = "", path) {
  Data <- lapply(X, function(y) {
    fread(file = file.path(path,
                           y))
  }) %>%
    # Bind list by position
    rbindlist(., use.names = F) %>%
    as.data.frame(.)
  #Column names refers to masses (averaged over all samples)
  if (grepl("\\+", Data$V1)[2] == T) {
    Y = "Mass_pos_"
  } else {
    Y = "Mass_neg_"
  }
  colnames(Data) <- Data[which(Data$V1 == "m/z"),
                         c(-1, -length(Data))] %>%
    apply(
      X = .,
      MARGIN = 2,
      FUN = function(X) {
        round(mean(X), digits = 3)
      }
    ) %>%
    c("", ., "Background") %>%
    paste0(Y, .)
  #Remove m/z rows
  Data <- Data[-which(Data[1] == "m/z"),]
  #Rownames sample - replicate id
  rownames(Data) <- Data[[1]] %>%
    gsub(".*@", "", .) %>%
    substr(x = .,
           start = 1,
           stop = nchar(.) - 13) %>%
    gsub("#", "", .) %>%
    gsub("May_|May-", "", .)
  Data[1] <- NULL
  #Normalization
  Backg <- Data[length(Data)]
  Data[length(Data)] <- NULL
  if (Normalization == "Sum") {
    #With sum of all peaks' intensities
    for (i in 1:length(Data[[1]])) {
      Data[i,] <- Data[i,] / mean(t(Data[i,]))
    }
  }
  if (Normalization == "Backg") {
    #with background intensity
    for (i in 1:length(Data[[1]])) {
      Data[i,] <- Data[i,] / Backg[[1]][[i]]
    }
  }
  return(Data)
}

## Data analysis ----
# Function to find out significant differences and masses
# which increased after treatment
YminusX <- function(X, Y) {
  Increased_Masses <- c()
  # Calculate the differences of the means and the mean values of 
  # the non-exposed samples
  MeanX <- apply(X, FUN = mean , MARGIN = 2)
  Diff <- apply(Y, FUN = mean, MARGIN = 2) - MeanX
  # Select the significant (>95%) differences and normalize them
  for (j in 1:length(Diff)) {
    if (t.test(Y[[j]], X[[j]])$p.value > 0.05) {
      Diff[[j]] <- NA
    }
    else {
      Diff[[j]] <- Diff[[j]] / MeanX[[j]]
      # Make a list of masses whose intensites increased significantly
      if (Diff[[j]] > 0) {
        Increased_Masses <- c(Increased_Masses, names(Diff)[[j]])
      }
    }
  }
  OUT <- list(Diff, Increased_Masses)
  names(OUT) <- c("Significant differences", "Increased masses")
  return(OUT)
}


## Random Forest analysis ----
RF <- function(X,
               seed = 123,
               mtry_step = 3,
               label = NULL) {
  # Create training and test data set
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
  
  # Top 50 models according to OOB error
  Param_Optimization <-
    hyper_grid[order(hyper_grid$OOB_error), ] %>% head(., 50)
  
  # Use best tuning parameters
  best_set <- hyper_grid[order(hyper_grid$OOB_error),][1,]
  
  # Re-run model with impurity-based variable importance
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
  
  # Re-run model with permutation-based variable importance
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
  # Most important variables according to impurtiy-based VI in distinguishing TPGs
  names(m_ranger_impurity$variable.importance) <- sub("_", "",
                                                      sub("Mass_", "",
                                                          names(
                                                            m_ranger_impurity$variable.importance
                                                          )))
  p1 <- vip::vip(m_ranger_impurity,
                 num_features = 20,
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
  
  # Most important variables according to permutation-based VI in distinguishing TPGs
  names(m_ranger_permutation$variable.importance) <- sub("_", "",
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

  # Access vip data
  var_imp <-
    vip::vip(m_ranger_impurity,
             num_features = 25,
             geom = NULL)
  RF_masses <- var_imp$data[[1]]
  
  # Predict for the test data
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
    plot_impurity = p1,
    plot_permutation = p2
  )
  return(O)
}

MDSplot_modified <- function (rf,
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
