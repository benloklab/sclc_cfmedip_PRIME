##################################################################
# SCLC Peripheral blood leukocyte methylation (PRIME) subtraction
# https://doi.org/10.1016/j.isci.2022.105487
#
# This script combines the CpG density data and the beta-values per windows
# and selects windows with high CpG density and low-PBL beta values
# Written by Sami Ul Haq, Nov 2021


# the following is the output from running step2.CpG.per.window.R
data.files <- Sys.glob("*num.cpgs.per.window.RData")

# holds all the cg per window data
num.cg.per.window <- c()

for(blah in data.files) {
  load(blah)
  num.cg.per.window <- c(num.cg.per.window, running.total.cgs.per.window)
}

# to save on memory
rm(running.total.cgs.per.window)

# genome wide CG density per window saved
save(num.cg.per.window, file="num.cg.per.window.RData")


# Num of CpG density per window is saved
summary(num.cg.per.window)
hist(num.cg.per.window, breaks=50)

# some output messages to get an idea of the number of windows
cat("\nMedian num of CG per 300bp window: ", median(num.cg.per.window), "\n")
# Num of CGs > 2 per 300bp window
cat("\nNum of CGs > 2 per 300bp window ", sum(num.cg.per.window > 2), "\n")
# Num of CGs > 5 per 300bp window
cat("\nNum of CGs > 5 per 300bp window ", sum(num.cg.per.window > 5), "\n")
# Num of CGs > 7 per 300bp window
cat("\nNum of CGs > 2 per 300bp window ", sum(num.cg.per.window > 7), "\n")
# Num of CGs > 9 per 300bp window
cat("\nNum of CGs > 10 per 300bp window ", sum(num.cg.per.window > 9), "\n")


#########################################################################################
# Using PBL depleted windows that had beta < 0.3 in SCLC PBLs
load("PBL.depleted.windows.beta_0.3.RData")

subsetted.windows <- num.cg.per.window[which(names(num.cg.per.window) %in% names(PBL.depleted.windows))]

# Num of CGs > 2 per 300bp window (~1.4M)
sum(subsetted.windows > 2)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 2)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.3.cg_2.RData")

# Num of CGs > 5 per 300bp window (~350k)
sum(subsetted.windows > 5)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 5)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.3.cg_5.RData")

# Num of CGs > 7 per 300bp window (~200k)
sum(subsetted.windows > 7)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 7)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.3.cg_7.RData")


# Num of CGs > 9 per 300bp window (~100k)
sum(subsetted.windows > 9)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 9)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.3.cg_9.RData")



#########################################################################################
# Using PBL depleted windows that had beta < 0.2 in SCLC PBLs
load("PBL.depleted.windows.beta_0.2.RData")

subsetted.windows <- num.cg.per.window[which(names(num.cg.per.window) %in% names(PBL.depleted.windows))]

# Num of CGs > 2 per 300bp window (~650k)
sum(subsetted.windows > 2)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 2)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.2.cg_2.RData")

# Num of CGs > 5 per 300bp window (~200k)
sum(subsetted.windows > 5)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 5)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.2.cg_5.RData")

# Num of CGs > 7 per 300bp window (~114k)
sum(subsetted.windows > 7)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 7)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.2.cg_7.RData")

# Num of CGs > 10 per 300bp window (~68k)
sum(subsetted.windows > 9)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 9)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.2.cg_9.RData")




#########################################################################################
# Using PBL depleted windows that had beta < 0.1 in SCLC PBLs
load("PBL.depleted.windows.beta_0.1.RData")

subsetted.windows <- num.cg.per.window[which(names(num.cg.per.window) %in% names(PBL.depleted.windows))]

# Num of CGs > 2 per 300bp window (~570k)
sum(subsetted.windows > 2)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 2)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.1.cg_2.RData")

# Num of CGs > 5 per 300bp window (~160k)
sum(subsetted.windows > 5)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 5)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.1.cg_5.RData")

# Num of CGs > 7 per 300bp window (~95k)
sum(subsetted.windows > 7)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 7)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.1.cg_7.RData")

# Num of CGs > 10 per 300bp window (~57k)
sum(subsetted.windows > 9)
PBL.depleted.windows <- subsetted.windows[which(subsetted.windows > 9)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.1.cg_9.RData")

