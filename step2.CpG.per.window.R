##################################################################
# SCLC Peripheral blood leukocyte methylation (PRIME) subtraction
# https://doi.org/10.1016/j.isci.2022.105487
#
# Step 2
# This script calculates the # of CGs (CpG density) per 300 bp window
# It takes an RData with file as command line argument
# Written by Sami Ul Haq

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

args <- commandArgs(trailingOnly=TRUE)

# RData file passed as 
the.file.to.use <- args[1]
load(the.file.to.use)

cat("\nworking with file", the.file.to.use, "\n\n")

# FYI, the RData file loads a vector called buncha.300bp.windows (change as needed)
buncha.windows <- buncha.300bp.windows

running.total.cgs.per.window <- c()

# The Number of CGs per window is calculated for all windows

debug.counter <- 0

cat("\nCalculating CGs per window\n")

for (each.window in buncha.windows) {
  # the window name is stored
  window.name <- each.window
  # the chromosome, start, and end is stored
  window.coordinates <- unlist(strsplit(window.name, split="[.]"))
  # a dataframe followed by a GRange is made (convenient data structure for subsequent analysis)
  temp.df <- data.frame(seqnames=window.coordinates[1],
                        start=window.coordinates[2],
                        end=window.coordinates[3])
  temp.grange <- makeGRangesFromDataFrame(temp.df)
  rm(temp.df)
  
  # the fragment of DNA associated with the above window is pulled
  frag.of.DNA <- getSeq(BSgenome.Hsapiens.UCSC.hg19, temp.grange, as.character=TRUE)
  
  num.of.cgs <- as.numeric(vcountPattern("CG",frag.of.DNA))
  running.total.cgs.per.window <- c(running.total.cgs.per.window, num.of.cgs)
  
  debug.counter <- debug.counter + 1
  
  if(debug.counter %% 1000 == 0) {
    cat("\nCalculating window ", debug.counter, "\n")
  }
}

cat("\nFinished calculating CGs per window\n")

names(running.total.cgs.per.window) <- buncha.windows
summary(running.total.cgs.per.window)

name.of.file.to.save <- paste0(the.file.to.use, ".num.cpgs.per.window.RData")

save(running.total.cgs.per.window, file=name.of.file.to.save)

