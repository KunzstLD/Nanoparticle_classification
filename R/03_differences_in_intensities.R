################################
#  DIFFERENCES IN INTENSITIES  #
################################

# Determine and list significant differences between exposed samples (marked with "F") and non-exposed samples
Diffs <- list()
Increased_Masses <- c()

for (i in c(1:10)){
  # Select the data corresponding to fulvic acid exposure and non-exposed
  Ful   <- Data[grep(rownames(Data), pattern = paste0("F",i,"-[1-9]")),]
  NoFul <- Data[grep(rownames(Data), pattern = paste0("^",i,"-[0-9]")),]
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
for (i in c(1:length(Diffs))){
  A <- data_frame(diff = Diffs[[i]], 
                  mass =  as.numeric(sub("_",".",sub("Mass_", "", names(Diffs[[i]])))),
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
write.csv(file = "Common masses", x = Common_masses)

# Cleanup
rm(Ful,NoFul,Diff,i,j,Means,A,Common_masses,Diff_plot)