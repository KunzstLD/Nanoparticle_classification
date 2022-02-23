# ___________________________________________________________________________
# DIFFERENCES IN INTENSITIES  ----
# Determine and list significant differences between exposed samples (marked with "F")
# and non-exposed samples
# ___________________________________________________________________________

# Read in nanoparticle masses
Masses_B <- readRDS(file.path(data_cache, "Masses_B.rds"))

# What has to be compared
Pairs <- data.frame(
  X = c("wafer", "^3", "^4", "9", "^3"),
  Y = c("^blank", "F3-", "F4-", "F9", "F3blank")
)

# Compare
Diffs <- list()
for (i in 1:nrow(Pairs)) {
  X <- Masses_B[grepl(Pairs$X[i], rownames(Masses_B)),]
  Y <- Masses_B[grepl(Pairs$Y[i], rownames(Masses_B)),]
  Diffs[[i]] <- YminusX(X, Y)
}
names(Diffs) <- c("blank-wafer","F3-3","F4-4","F9-9","F3blank-3")

# Export lists of fragments
for (i in 2:4) {
  write.csv(Diffs[[i]]$'Significant differences',
            file = file.path(data_out, paste0(
              names(Diffs)[i], "_Significant_differences.csv"
            )))
  write.csv(Diffs[[i]]$'Increased Masses_B',
            file = file.path(data_out, paste0(names(Diffs)[i], "_Increased_Masses_B.csv")))
}

# Preparing for plotting
Diff_plot <- data.frame()
for (i in c(1:length(Diffs))){
  A <- tibble(diff = Diffs[[i]]$`Significant differences` %>% as.numeric(.),
              mass = Diffs[[i]]$`Significant differences` %>% 
                     names(.) %>% 
                     sub("Mass_neg_|Mass_pos_","",.) %>% 
                     as.numeric(.),
              polarity = c(
                            Diffs[[i]]$`Significant differences` %>% names(.) %>%
                              grep("neg",.) %>% length(.) %>%
                              rep("neg",.),
                            Diffs[[i]]$`Significant differences` %>% names(.) %>%
                              grep("pos",.) %>% length(.) %>%
                              rep("pos",.)
                           ),
              sample = rep(names(Diffs)[[i]], length(Diffs[[i]]$`Significant differences`)) %>% 
                       as.factor(.)
  )
  Diff_plot <- rbind(Diff_plot, A)
}

# Plot
ggplot(Diff_plot, aes(x = mass, y = diff, color = polarity)) +
  geom_point(size = 0.8, shape = 21) +
  facet_wrap( ~ sample,
              scale = "free_y",
              ncol = 2) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = element_blank(),
    x = "m/z",
    y = "Normalized intensity differences",
    linetype = NULL
  )

# Save R object for differences 
saveRDS(object = Diffs, file = file.path(data_cache, "Diffs.rds"))
