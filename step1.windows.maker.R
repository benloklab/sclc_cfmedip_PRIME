##################################################################
# SCLC Peripheral blood leukocyte methylation (PRIME) subtraction
# https://doi.org/10.1016/j.isci.2022.105487
#
# Step 1
# This script creates genome-wide 300bp sized windows spanning chromosome 1-22
# Written by Sami Ul Haq 2020

library(rtracklayer)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

chromosomes.selected <- c(1:22)

# creates 300bp whole genome (chr 1-22)
genome_300bp <- data.frame(genomeBlocks(BSgenome.Hsapiens.UCSC.hg19, chrs=seqnames(BSgenome.Hsapiens.UCSC.hg19)[chromosomes.selected], width=300))
genome_300bp$name_format <- paste(genome_300bp$seqnames, genome_300bp$start, genome_300bp$end, sep=".")
genome_300bp <- makeGRangesFromDataFrame(genome_300bp, keep.extra.columns = TRUE)

buncha.300bp.windows <- genome_300bp$name_format

save(buncha.300bp.windows, file="genome.wide.300bp.windows.RData")