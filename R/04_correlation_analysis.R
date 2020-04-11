############################
##  CORRELATION ANALYSIS  ##
############################

# Remove S2 (no increased masses)
Diffs_pos$S2 <- NULL

# Build a list of network plots
list_plot <- list()
for (i in c(1:length(Diffs_pos))){
  # Determine the correlation matrix for the increased masses of the sample in the whole data set
  Cor.matrix <- cor(Data[,names(Diffs_pos[[i]])], method = "pearson")
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

#Clean up
rm(i,Cor.matrix,list_plot)