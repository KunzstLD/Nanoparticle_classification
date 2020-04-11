# Import data
# List all files to be read in
files <- list.files(path = file.path(data_in),
                    pattern = ".*TXT")

# Read in
Data <- lapply(files, function(y) {
  fread(file = file.path(data_in,
                         y),
        drop = c("m/z", "V200"))
})

# Bind list by position
Data <- rbindlist(Data, use.names = F)

# Change colnames
names(Data)[52] <- "52,01"
setnames(Data,
         old = names(Data),
         new = paste0("Mass_", gsub("\\,", "_" , names(Data))))

# all integer columns converted to numeric
cols <- names(Filter(is.integer, Data))
Data[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

# Add id as rownames
# (data.table does not support rownames, hence we convert Data back to a data.frame)
setDF(Data)
IDs <- paste0(rep(1:10, each = 9), "-", c(1:9))
for (i  in 1:10) {
  IDs <- c(IDs, paste0("F", i, "-", c(1:9)))
}
rownames(Data) <- IDs

# Clean up
rm(IDs, i, files)
