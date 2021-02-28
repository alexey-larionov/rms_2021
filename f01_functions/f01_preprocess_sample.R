# Function to preprocess sample in RMS progect

# Alexey Larionov, 25 Feb 2021

# Summary
#--------

# Preprocesses sample.gr:  
# - Removes variants on non-standard cheromosomes
# - Keeps only variants overlapping with the genes
# - Removes common variants (gnomAD MAF > 0.05)

# Thakes two GRanges
# Outputs one GRange

# Assumes that sample.gr and genes.gr have the same Seqinfo
# Assumes that sample.gr contains gnomAD AF-s

# Required libraries
#-------------------

require(GenomicRanges)

# Function
#---------

preprocess_sample.rms <- function(sample.gr, genes.gr){

  # Keep only variants on standard chromosomes
  sample.gr <- keepStandardChromosomes(sample.gr, pruning.mode="coarse")

  # Overlap with the genes list
  sample.gr <- subsetByOverlaps(sample.gr, genes.gr, type="any")

  # Exclude common variants (MAF > 0.05)
  sample.gr <- sample.gr[sample.gr$gnomAD_genome_ALL <= 0.05 |
                         sample.gr$gnomAD_genome_ALL >= 0.95 |
                         is.na(sample.gr$gnomAD_genome_ALL)]
  
  # Return result
  return(sample.gr)

}

# Pregress report
cat("Loaded function preprocess_sample.rms()","\n")
