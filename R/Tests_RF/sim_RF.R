# TODO: How many samples for the individual sun screen extracts?
# Idea: Make up training dataset (Use bootstrap samples from our data?)
# -> Stepwise: 20 (original), 30, 40, 50, 60, 70, 80, 90, 100, ...
# Evaluate Model skill: Prediction accuracy?
# Look at VI results
# -> Next 2-3 weeks

# ____________________________________________________________
# What's the sample size for stabilizing the prediction error?
# What's the sample size for stabilizing permutation/impurtiy-based VI measures?
# ____________________________________________________________

# libraries
# data cleaning/ data transformation
library(data.table) 
library(dplyr)
library(tidyverse)

# random forest
library(rsample)
library(ranger)
library(randomForest)

# data: load Data
# load 01_preprocessing script
# Add column indicating presence of fulvic acid
Data_RF <-
  cbind(Data,
        fulvic_acid = as.factor(ifelse(grepl("F", rownames(
          Data
        )),
        1, 0)))

# Remove sample 2 (no sorption)
Data_RF <- Data_RF[!grepl("F2|2", rownames(Data)),]

# save rownames
samples <- rownames(Data_RF)

# create dt with sample and smaple_id columns
setDT(Data_RF)
Data_RF[, nanopart := samples]
Data_RF[, nanopart_ids := sub("\\-.*", "", nanopart)]

# RF function 
fun_optimized <- function(X,
                #seed = 123,
                mtry_step = 1,
                label = NULL) {
  # create training and test data set
  # set.seed(seed)
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
  # m_ranger_permutation <- ranger(
  #   formula         = fulvic_acid ~ .,
  #   data            = X_train,
  #   num.trees       = n_features * 10,
  #   mtry            = best_set$mtry,
  #   min.node.size   = best_set$node_size,
  #   sample.fraction = best_set$sample_size,
  #   importance      = 'permutation',
  #   verbose         = FALSE,
  #   seed            = 123
  # )
  # most important variables according to impurtiy-based VI in distinguishing masses
  names(m_ranger_impurity$variable.importance) <- sub("_", ".",
                                                      sub("Mass_", "",
                                                          names(
                                                            m_ranger_impurity$variable.importance
                                                          )))
  
  #### Predict for the test data ####
  pred_class <-
    predict(m_ranger_impurity, X_test[, -c("fulvic_acid")])
  
  # Assess performance on test data
  Confusion_matrix <-
    caret::confusionMatrix(factor(pred_class$predictions),
                           factor(X_test$fulvic_acid))
  
  output <- list(m_ranger_impurity$prediction.error, 
                 pred_class$num.samples, 
                 Confusion_matrix$overall[["Accuracy"]])
  names(output) <- c("prediction.error", "num.samples.test", "accuracy.test")
  return(output)
}

fun_base <- function(X,
                     #seed = 123,
                     mtry_step = 1,
                     label = NULL) {
  # create training and test data set
  # set.seed(seed)
  split <- initial_split(X, prop = .7, strata = "fulvic_acid")
  X_train <- training(split)
  X_test  <- testing(split)
  
  # number of features
  n_features <- length(setdiff(names(X_train), "fulvic_acid"))
  
  # Grid for different parameters
  
  # execute rf with default
  m_ranger_impurity <- ranger(
    formula         = fulvic_acid ~ .,
    data            = X_train,
    num.trees       = n_features * 10,
    mtry            = 2,
    min.node.size   = 1,
    sample.fraction = 0.63,
    importance      = 'impurity',
    verbose         = FALSE,
    seed            = 123
  )
  
  # most important variables according to impurtiy-based VI in distinguishing masses
  names(m_ranger_impurity$variable.importance) <- sub("_", ".",
                                                      sub("Mass_", "",
                                                          names(
                                                            m_ranger_impurity$variable.importance
                                                          )))
  
  #### Predict for the test data ####
  pred_class <-
    predict(m_ranger_impurity, X_test[, -c("fulvic_acid")])
  
  # Assess performance on test data
  Confusion_matrix <-
    caret::confusionMatrix(factor(pred_class$predictions),
                           factor(X_test$fulvic_acid))
  
  output <- list(
    m_ranger_impurity$prediction.error,
    pred_class$num.samples,
    Confusion_matrix$overall[["Accuracy"]]
  )
  names(output) <-
    c("prediction.error", "num.samples.test", "accuracy.test")
  return(output)
}

# ____________________________________________________________________
##### 1) Create bootstrap samples from S4 sample (40 predictors) ####
# Does always produce the same results and 
# RF overfits 
# Probably due to creation of samples
# ____________________________________________________________________

# 
## load increased masses
# Diffs_pos <- readRDS(file = file.path(data_out, "Diffs_pos.rds"))
# Diffs_pos$S2 <- NULL

# two questions: 
# compare to test data -> and then evaluate prediction error
# USE the Relative root mean squared error (RRE)

# sample_sizes <- seq(20,100, by = 10)
# result <- list()
# 
# for(n in sample_sizes) {
#   # sample training data
#   set.seed(1234)
#   
#   # create subsets and select only masses per sample
#   # that showed increasing masses
#   cols <- c(names(Diffs_pos[[3]]), "fulvic_acid")
#   nofulvic <- Data_RF[nanopart_ids == "4",
#                       lapply(.SD, function(y)
#                         sample(y, size = n, replace = TRUE)),
#                       .SDcols = cols]
#   
#   fulvic <- Data_RF[nanopart_ids == "F4",
#                     lapply(.SD, function(y)
#                       sample(y, size = n, replace = TRUE)),
#                     .SDcols = cols]
#   # bind together
#   nano_data <- rbind(nofulvic,
#                      fulvic)
#   result[[paste(n)]] <- fun(X = nano_data)
# }

# system.time(fun(X = nano_data), 
#             gcFirst = FALSE)

# __________________________________________________________________________
#### 2) Pool all masses irrespective of nanoparticle origin and sample #####
# __________________________________________________________________________
# Why does replicate always give the same results?
# obj <- replicate(5, {
#   test <- fun(X = nano_data)
#   test
# })

# Remove masses whose intensities decreased
# Increased_masses <- readRDS(file = file.path(data_out, "Increased_Masses.rds"))
Increased_masses <- readRDS(file = file.path(getwd(), "Cache", "Increased_Masses.rds"))
Data_RF_inc <- copy(Data_RF[, .SD, .SDcols = c(Increased_masses, 
                                               "fulvic_acid", 
                                               "nanopart",
                                               "nanopart_ids")])

# sample from each mass 
cols <- grep("Mass.*|fulvic.*", names(Data_RF_inc), value = TRUE)


# Trial 1: Six datasets, 5 replicates, 1 time RF 
# output <- list()
# sample_sizes <- seq(5, 10, by = 5)
# 
# for(n in sample_sizes) {
#   dat <- replicate(5, {
#     nofulvic <-
#       Data_RF_inc[nanopart_ids %in% c("1":"10"), lapply(.SD, function(y)
#         sample(y, size = n, replace = TRUE)), .SDcols = cols]
# 
#     fulvic <-
#       Data_RF_inc[nanopart_ids %like% "F.*", lapply(.SD, function(y)
#         sample(y, size = n, replace = TRUE)),
#         .SDcols = cols]
# 
#     nano_data <- rbind(nofulvic,
#                        fulvic)
#     nano_data
#   },
#   simplify = FALSE)
# 
#   output[[paste(n)]] <- dat
# }
#  
# nano_data_bt <-
#   lapply(output, function(y)
#     rbindlist(y, idcol = "id")) %>%
#   rbindlist(., idcol = "file") %>%
#   .[, true_id := paste0(file, "_", id)] %>%
#   .[, c("id", "file") := NULL]

# first approach using classical lapply
# system.time(result_rf_sim <- lapply(output, function(nested) {
#   lapply(nested, function(y) {
#     replicate(2, {
#       fun(y)
#     })
#   })
# }), 
# gcFirst = FALSE)

# first test run: 1 dataset, 5 replicates, 1 RF
# took ~ 3 min
# test <- rbindlist(output[[1]], idcol = "id")
# system.time(test_res <- test[, fun(X = .SD),
#                              .SDcols = cols,
#                              by = "id"])

# second run: 6 datasets, 5 replicates, 1 RF
# should take ~ 18 min
# Reality: took ~ 29 mins
# base::print(system.time(res_nano_6 <- nano_data_bt[, fun(X = .SD),
#                              .SDcols = cols,
#                              by = "true_id"]))
# print(res_nano_6)


# Trial 2: 12 datasets, 5 replicates, 1 time RF 
# output <- list()
# sample_sizes <- seq(5, 60, by = 5)
# 
# for(n in sample_sizes) {
#   dat <- replicate(5, {
#     nofulvic <-
#       Data_RF_inc[nanopart_ids %in% c("1":"10"), lapply(.SD, function(y)
#         sample(y, size = n, replace = TRUE)), .SDcols = cols]
#     
#     fulvic <-
#       Data_RF_inc[nanopart_ids %like% "F.*", lapply(.SD, function(y)
#         sample(y, size = n, replace = TRUE)),
#         .SDcols = cols]
#     
#     nano_data <- rbind(nofulvic,
#                        fulvic)
#     nano_data
#   },
#   simplify = FALSE)
#   
#   output[[paste(n)]] <- dat
# }
# 
# nano_data_bt_12 <-
#   lapply(output, function(y)
#     rbindlist(y, idcol = "id")) %>%
#   rbindlist(., idcol = "file") %>%
#   .[, true_id := paste0(file, "_", id)] %>%
#   .[, c("id", "file") := NULL]
# 
# 
# # third run: 12 datasets, 5 replicates, 1 RF
# # should take ~ 36 min
# base::print(system.time(res_nano_12 <- nano_data_bt_12[, fun(X = .SD),
#                                                    .SDcols = cols,
#                                                    by = "true_id"]))
# print(res_nano_12)


# Trial 3: 12 datasets, 10 replicates, 1 time RF, RF base ----
output <- list()
sample_sizes <- seq(5, 60, by = 5)
set.seed(123)

for(n in sample_sizes) {
  dat <- replicate(10, {
    nofulvic <-
      Data_RF_inc[nanopart_ids %in% c("1":"10"), lapply(.SD, function(y)
        sample(y, size = n, replace = TRUE)), .SDcols = cols]
    
    fulvic <-
      Data_RF_inc[nanopart_ids %like% "F.*", lapply(.SD, function(y)
        sample(y, size = n, replace = TRUE)),
        .SDcols = cols]
    
    nano_data <- rbind(nofulvic,
                       fulvic)
    nano_data
  },
  simplify = FALSE)
  
  output[[paste(n)]] <- dat
}

nano_data_bt_12_10 <-
  lapply(output, function(y)
    rbindlist(y, idcol = "id")) %>%
  rbindlist(., idcol = "file") %>%
  .[, true_id := paste0(file, "_", id)] %>%
  .[, c("id", "file") := NULL]

# third run: 12 datasets, 10 replicates, 1 RF
base::print(system.time(res_nano_12_10 <-
                          nano_data_bt_12_10[, fun_base(X = .SD),
                                            .SDcols = cols,
                                            by = "true_id"]))

#### Simulation with base RF ============================================
RF_base_sim <- replicate(100, nano_data_bt_12_10[, fun_base(X = .SD),
                                                 .SDcols = cols,
                                                 by = "true_id"],
                         simplify = FALSE)

#### Analysis ============================================================
RF_base_sim <- rbindlist(RF_base_sim, idcol = "file")
setnames(RF_base_sim, 
         "file",
         "RF_run")

# create col with sample sizes
RF_base_sim[, sample.size := as.numeric(sub("\\_.*", "", true_id))]

# plot
melt(RF_base_sim,
     measure.vars = c("prediction.error", "accuracy.test")) %>%
  ggplot(., aes(
    x = as.factor(sample.size),
    y = value,
    color = variable
  )) +
  geom_boxplot(size = 0.6, alpha = 0.5) +
  geom_hline(yintercept = 0.05,
             color = "grey",
             linetype = "dashed") +
  #scale_color_hue(l=60, c=45)+
  scale_color_manual(values = c("#56B4E9", "#E69F00")) + # "#56B4E9"
  theme_classic() +
  labs(x = "Sample size", y = NULL) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(family = "Roboto Mono", size = 11),
    axis.text.y = element_text(family = "Roboto Mono", size = 11),
    panel.grid = element_blank()
  ) +
  annotate(
    "text",
    x = 11.5,
    y = 0.12,
    family = "Poppins",
    size = 3,
    color = "gray20",
    label = "prediction error or \n accuracy of 0.05"
  )
ggsave(
  file.path(data_out, "RF_sim_results.png"),
  width = 25,
  height = 12,
  units = "cm"
)





