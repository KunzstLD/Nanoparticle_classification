# ___________________________________________________________________________
# Import data ----
# The first rows containing "null" and "Area statistics" have to be deleted.
# ___________________________________________________________________________

# List all files to be read in
files_neg <- list.files(path = file.path(data_in),
                        pattern = "neg.*TXT")
files_pos <- list.files(path = file.path(data_in),
                        pattern = "pos.*TXT")

# Combine positive and negative polarity
Masses <-
  cbind(as.data.frame(Import(files_neg, Normalization = "Sum")),
        as.data.frame(Import(files_pos, Normalization = "Sum")))

# Save as rds object to Cache folder
saveRDS(object = Masses, file = file.path(data_cache, "Masses.rds"))
