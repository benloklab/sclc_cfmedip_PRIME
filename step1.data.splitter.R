##################################################################
# SCLC Peripheral blood leukocyte methylation (PRIME) subtraction
# https://doi.org/10.1016/j.isci.2022.105487
#
# Step 1
# This script splits the PBL depleted windows data 
# file into smaller chunks, as this makes computation easier
# Written by Sami Ul Haq Oct 2021

# Data needs to be downloaded from the following Zenodo repository:
# https://zenodo.org/record/8286762
# DOI: 10.5281/zenodo.8286762
#
# Full citation:
# Ul Haq, Sami, Tsao, Ming S., Cabanero, Michael, de Carvalho, Daniel, Liu, 
# Geoffrey, Bratman, Scott V., & Lok, Benjamin H. (2022). 
# Processed counts data of cfMeDIP-seq profiles of small cell 
# lung cancer patients [Data set]. In Cell Press iScience (Vol. 25, Number 12). 
# Zenodo. https://doi.org/10.5281/zenodo.8286762


# the "medestrand.300bp.matrix.of.pbls.RData" file can be found in the repository described in the above comments
load("medestrand.300bp.matrix.of.pbls.RData")
PBL.depleted.windows <- rownames(matrix.of.pbls)

# The data file will be used to calculate CpG density. To do this efficiently, the data is split up in smaller parts
# In the code below, the code is divided into 2million windows
subsetted <- PBL.depleted.windows[1:2000000]
save(subsetted, file="PBL.depleted.windows.part1.RData")

subsetted <- PBL.depleted.windows[2000001:4000000]
save(subsetted, file="PBL.depleted.windows.part2.RData")

subsetted <- PBL.depleted.windows[4000001:6000000]
save(subsetted, file="PBL.depleted.windows.part3.RData")

subsetted <- PBL.depleted.windows[6000001:8000000]
save(subsetted, file="PBL.depleted.windows.part4.RData")

subsetted <- PBL.depleted.windows[8000001:length(PBL.depleted.windows)]
save(subsetted, file="PBL.depleted.windows.part5.RData")

