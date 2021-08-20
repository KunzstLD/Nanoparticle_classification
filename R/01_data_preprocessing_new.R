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

#Removing trends
#Correct <- function(X, names, numb_rep, df = 3){
 # for (i in names){
   # A <- X[which(grepl(paste0("^",i),rownames(X))==T),]
   # for (j in 1:length(A)){
     # Trend <- A[[j]] %>% smooth.spline(x=c(1:numb_rep),y=.,df=df)
     # A[[j]] <- median(A[[j]])+A[[j]]-predict(Trend,x=c(1:numb_rep))$y
   # }
   # X[which(grepl(paste0("^",i),rownames(X))==T),] <- A
 # }
 # return(X)
#}
#Data <- Correct(Data, 
                #names = c("wafer-","3-", "4-", "9-", "blank-", "F3-", "F3blank-", "F4-", "F9-"), 
               # numb_rep = 36, 
               # df = 3)

# Optimize correction
Sample_names <-
  c("wafer-",
    "3-",
    "4-",
    "9-",
    "blank-",
    "F3-",
    "F3blank-",
    "F4-",
    "F9-")
for (i in Sample_names) {
  A <- Masses[which(grepl(paste0("^", i), rownames(Masses)) == T), ]
  par(
    mfrow = c(3, 4),
    mar = c(4, 1, 3, 2),
    oma = c(1.5, 2, 2, 1)
  )
  for (j in c(1:4, 201:204, 501:504)) {
    Trend <- smooth.spline(x = c(1:36),
                           y = A[[j]],
                           df = 3)
    plot(A[[j]],
         main = colnames(Masses)[j],
         ylab = "",
         xlab = "")
    lines(Trend)
    points(median(A[[j]]) + A[[j]] - predict(Trend, x = c(1:36))$y,
           col = "red",
           pch = 19)
  }
  par(
    fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0),
    new = TRUE
  )
  plot(
    0,
    0,
    type = 'l',
    bty = 'n',
    xaxt = 'n',
    yaxt = 'n'
  )
  legend(
    "bottom",
    legend = c("uncorrected", "corrected"),
    col = c("black", "red"),
    horiz = T,
    lwd = 5,
    xpd = T,
    cex = 1,
    seg.len = 1
  )
  mtext(
    sub("-", "", i),
    cex = 1,
    side = 3,
    line = -2,
    font = 2
  )
}

# Save as rds object to Cache folder
saveRDS(object = Masses, file = file.path(data_cache, "Masses.rds"))
