# ___________________________________________________________________________
#  RANDOM FOREST  ----
# ___________________________________________________________________________

# Read in nanoparticle data
Masses_B <- readRDS(file.path(data_cache, "Masses_B.rds"))
Diffs <- readRDS(file.path(data_cache, "Diffs.rds"))


Samples_names <- c("3-","4-","9-")
Results_RF_Extracts <- list()
for (Sample in ) {
  # Select samples to be processed
  Inc_Masses_B <-
    Diffs[[grep(Sample, names(Diffs))]]$`Increased masses`
  Data_RF <- Masses_B[grepl(Sample , rownames(Masses_B)), Inc_Masses_B]
  
  # Add column indicating presence of fulvic acid
  Data_RF$fulvic_acid <-
    as.factor(ifelse(grepl("F", rownames(Data_RF)),
                     1, 0))
  
  # Carry out Random forest
  Results_RF_Extracts[[Sample]] <-
    RF(Data_RF, label = paste0("S", sub("-\\d+", "", rownames(Data_RF)[1])))
  
  # Plot importance plots
  p1 <- Results_RF_Extracts[[Sample]]$plot_impurity
  p2 <- Results_RF_Extracts[[Sample]]$plot_permutation
  png(file = file.path(data_out, paste0(
    "RF_importances_", "S", sub("\\-", "", Sample), ".png"
  )))
  gridExtra::grid.arrange(p1, p2, nrow = 1)
  dev.off()
  
  # Plot Proximity map
  model1 <-
    randomForest(
      fulvic_acid ~ .,
      data = Data_RF,
      ntree = length(Data_RF)*10,
      mtry = Results_RF_Extracts[[Sample]]$Parameters_optimization[1,1],
      nodesize = Results_RF_Extracts[[Sample]]$Parameters_optimization[1,2],
      importance = T,
      proximity = T
    )

  # save proximity map
  pdf(file = file.path(data_out, paste0("RF_proxmap_", "S", sub("\\-", "", Sample),".pdf")))
  MDSplot_modified(model1,
                   Data_RF$fulvic_acid,
                   palette = c("blue", "red"),
                   k = 2)
  dev.off()
}
saveRDS(object = Results_RF_Extracts, file = file.path(data_out, "RF_results.rds"))