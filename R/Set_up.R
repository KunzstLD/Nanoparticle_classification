# libraries
# data cleaning/ data transformation
library(data.table) 
library(dplyr)
library(tidyverse)

# Decision tree
library(dtree)
library(rpart.plot)
library(rpart.utils)

# random forest
library(rsample)
library(ranger)
library(randomForest)

# correlation analysis
library(corrr)

# clustering
library(optpart)

# plotting
library(ggplot2)
library(cowplot)
library(grid)
library(ggdendro)
library(plotly)

# specify input data path
# needs to be set by the user
data_in <- "./Raw_data"

data_out <- "./Output"

data_scr <- "./R"

data_cache <- "./Cache"

# Load functions script
source(file.path(data_scr, "functions.R"))
