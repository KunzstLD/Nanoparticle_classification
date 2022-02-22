# ___________________________________________________________________________
# Import data ----
# The first rows containing "null" and "Area statistics" have to be deleted.
# ___________________________________________________________________________

## First dataset ----
# Part A

# List all files to be read in
files <- list.files(path = file.path(data_in, "A"), pattern = ".*TXT")

# Read in
Data <- lapply(files, function(y) {
  fread(file = file.path(data_in, "A", y),
        drop = c("m/z", "V200"))
})

# Bind list by position
Data <- rbindlist(Data, use.names = F) 

# Change colnames
#colnames(Data)[52] <- "52,01"
setnames(Data,
         old = names(Data),
         new = paste0("Mass_", gsub("\\,", "_" , names(Data))))

# Convert to numeric
Data <- as.data.frame(apply(Data, 
                            MARGIN = c(1,2), 
                            FUN = function(x){as.numeric(x)}))     

# Add id column
IDs <- c()
for (i  in 1:10){IDs <- c(IDs, paste0(i,"-", c(1:9)))}
for (i  in 1:10){IDs <- c(IDs, paste0("F",i,"-", c(1:9)))}
rownames(Data) <- IDs
saveRDS(object = Data, file = file.path(data_cache, "Data.rds"))


## Second dataset ---- 
# Part B
# (3 sunscreens)
# List all files to be read in
files_neg <- list.files(path = file.path(data_in, "B"),
                        pattern = "neg.*TXT")
files_pos <- list.files(path = file.path(data_in, "B"),
                        pattern = "pos.*TXT")

# Combine positive and negative polarity
Masses <-
  cbind(as.data.frame(Import(files_neg, Normalization = "Sum", path = file.path(data_in, "B"))),
        as.data.frame(Import(files_pos, Normalization = "Sum", path = file.path(data_in, "B"))))

# Save as rds object to Cache folder
saveRDS(object = Masses, file = file.path(data_cache, "Masses.rds"))
