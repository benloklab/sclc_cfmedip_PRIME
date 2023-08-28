##################################################################
# SCLC Peripheral blood leukocyte methylation (PRIME) subtraction
# https://doi.org/10.1016/j.isci.2022.105487
#
# This script takes a buncha MeDIP-library BAM files of PBLs and
# calculates a beta-value for each 300bp window (converts counts data to beta values)
# Written By Sami Ul Haq, Sept-14-2021

library("MeDEStrand")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BiocGenerics")
library("S4Vectors")
library("IRanges")
library("GenomicRanges")
library("Biostrings")
library("XVector")
library("rtracklayer")
library("Rsamtools")
library("MEDIPS")
library("Repitools")

# parameters for binning genome
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq = 1
extend = 200
shift = 0
ws = 300
chr.select = paste0('chr', c(1:22) )

args <- commandArgs(trailingOnly=TRUE)
directory.to.use <- args[1]

bam.files.names <- Sys.glob(paste0(directory.to.use, "/*bam"))

# creates 300bp tiled genome for chromosomes 1 - 22
genome_300bp <- data.frame(genomeBlocks(BSgenome.Hsapiens.UCSC.hg19, chrs=seqnames(BSgenome.Hsapiens.UCSC.hg19)[1:22], width=300))
genome_300bp$name_format <- paste(genome_300bp$seqnames, genome_300bp$start, genome_300bp$end, sep=".")

# this matrix will hold all medestrand beta values
matrix.of.pbls <- matrix(nrow=length(genome_300bp$name_format))
# the rownames correspond to windows
rownames(matrix.of.pbls) <- genome_300bp$name_format

# each bam file is examined and the beta-value is infered using MeDEStrand and then added as a column to the holder matrix
for(each.bam.file in bam.files.names) {
  MeDIP_seq = MeDEStrand.createSet(file=each.bam.file, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr.select, paired = T)
  
  #  count CpG pattern in the bins
  CS = MeDEStrand.countCG(pattern="CG", refObj=MeDIP_seq)
  result.methylation = MeDEStrand.binMethyl(MSetInput = MeDIP_seq, CSet = CS, Granges = FALSE)
  
  # each sample is added
  matrix.of.pbls <- cbind(matrix.of.pbls, result.methylation)
}

# the first column of the matrix contains NAs
matrix.of.pbls <- matrix.of.pbls[,-c(1)]

# the colnames are the same as bam file names 
colnames(matrix.of.pbls) <- bam.files.names


save(matrix.of.pbls, file="medestrand.300bp.matrix.of.pbls.RData")


