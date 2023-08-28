##################################################################
# SCLC Peripheral blood leukocyte methylation (PRIME) subtraction
# https://doi.org/10.1016/j.isci.2022.105487
# 
# This script finds hypo-methylated windows in PBL samples (PBL-depleted windows)
# Written by Sami Ul Haq, Sept 15, 2021
#

# loads the matrix containing MeDEStrand converted beta-values for the PBLs
# this RData file loads a matrix called "matrix.of.pbls"
load("medestrand.300bp.matrix.of.pbls.RData")

# this examines the median beta values per window
median.beta.vals <- apply(matrix.of.pbls, MARGIN=1, FUN = median)
names(median.beta.vals) <- rownames(matrix.of.pbls)

# Filter out blacklist windows (V2 ENCODE-blacklist windows)
load("hg19.encode.blacklist.windows.RData")
# matrix is removed for blacklist windows
median.beta.vals <- median.beta.vals[ !(names(median.beta.vals) %in% hg19.encode.blacklist.windows) ]

# selects windows with beta values less than 0.3
PBL.depleted.windows <- median.beta.vals[which(median.beta.vals < 0.3)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.3.RData")

# selects windows with beta values less than 0.2
PBL.depleted.windows <- median.beta.vals[which(median.beta.vals < 0.2)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.2.RData")

# selects windows with beta values less than 0.1
PBL.depleted.windows <- median.beta.vals[which(median.beta.vals < 0.1)]
save(PBL.depleted.windows, file="PBL.depleted.windows.beta_0.1.RData")
