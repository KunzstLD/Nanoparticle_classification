# ___________________________________________________________________________
# Import data ----
# ___________________________________________________________________________

## First dataset ----
# Part A
# Analysis on sorption patterns

# List all files to be read in
files <- list.files(
  path = file.path(data_in, "A_SorptionPatterns_Sunsreens_1_10"),
  pattern = ".*TXT"
)

# Read in
Masses_A <- lapply(files, function(y) {
  fread(file = file.path(data_in, "A_SorptionPatterns_Sunsreens_1_10", y),
        drop = c("m/z", "V200"))
})

# Bind list by position
Masses_A <- rbindlist(Masses_A, use.names = F) 

# Change colnames
# colnames(Data)[52] <- "52,01"
setnames(Masses_A,
  old = names(Masses_A),
  new = paste0("Mass_", gsub("\\,", "_", names(Masses_A)))
)

# Convert to numeric
Masses_A <- as.data.frame(apply(Masses_A,
  MARGIN = c(1, 2),
  FUN = function(x) {
    as.numeric(x)
  }
))

# Add id column
IDs <- c()
for (i  in 1:10){IDs <- c(IDs, paste0(i,"-", c(1:9)))}
for (i  in 1:10){IDs <- c(IDs, paste0("F",i,"-", c(1:9)))}
rownames(Masses_A) <- IDs
saveRDS(
  object = Masses_A,
  file = file.path(data_cache, "Masses_A.rds")
)

## Second dataset ---- 
# Analysis on most important masses for sorption
# Part B
# (3 sunscreens)
# List all files to be read in
files_neg <- list.files(path = file.path(data_in, "B_ImportantMasses_Sunscreens_3_4_9"),
                        pattern = "neg.*TXT")
files_pos <- list.files(path = file.path(data_in, "B_ImportantMasses_Sunscreens_3_4_9"),
                        pattern = "pos.*TXT")

# Combine positive and negative polarity
Masses_B <-
  cbind(
    as.data.frame(Import(files_neg, Normalization = "Sum", path = file.path(data_in, "B_ImportantMasses_Sunscreens_3_4_9"))),
    as.data.frame(Import(files_pos, Normalization = "Sum", path = file.path(data_in, "B_ImportantMasses_Sunscreens_3_4_9")))
  )

# Save as rds object to Cache folder
saveRDS(object = Masses_B, file = file.path(data_cache, "Masses_B.rds"))
