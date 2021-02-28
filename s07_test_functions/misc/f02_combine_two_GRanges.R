# Function to combine two GRanges in RMS progect

# Alexey Larionov, 25 Feb 2021

# Summary
#--------

# Takes two GRanges and a Data-Frame.

# The Data-Frame contains samples descriptions:  
# Sample data for gr1 must be in record1 and 
# sample data for gr2 must be in record2 in the dataframe.  

# Outputs SummarizedExperiment with data merged from the GRanges.


# Required libraries
#-------------------

require(GenomicRanges)
require(SummarizedExperiment)

# Function
#---------

combine_two_GRanges.rms <- function(gr1, gr2, samples.df){
  
  # Add key to gr1
  #---------------
  
  # Select data for the key
  gr1_chr <- as.character(gr1@seqnames)
  gr1_start <- as.data.frame(gr1@ranges)$start
  gr1_end <- as.data.frame(gr1@ranges)$end
  gr1_ref <- gr1$Ref
  gr1_alt <- gr1$Alt
  
  # Make and add key to gr1
  gr1_var_id <- paste(gr1_chr, gr1_start, gr1_end, gr1_ref, gr1_alt, sep="_")
  gr1$var_id <- gr1_var_id
  
  # Add key to GRanges names
  names(gr1) <- gr1_var_id
  
  # Add key to gr2
  #---------------
  
  # Select data for the key
  gr2_chr <- as.character(gr2@seqnames)
  gr2_start <- as.data.frame(gr2@ranges)$start
  gr2_end <- as.data.frame(gr2@ranges)$end
  gr2_ref <- gr2$Ref
  gr2_alt <- gr2$Alt
  
  # Make and add key to gr2
  gr2_var_id <- paste(gr2_chr, gr2_start, gr2_end, gr2_ref, gr2_alt, sep="_")
  gr2$var_id <- gr2_var_id
  
  # Add key to GRanges names
  names(gr2) <- gr2_var_id
  
  # Select subsets
  #---------------
  
  unique_gr1 <- setdiff(gr1$var_id, gr2$var_id)
  unique_gr2 <- setdiff(gr2$var_id, gr1$var_id)
  intersect_gr1_gr2 <- intersect(gr1$var_id, gr2$var_id)
  union_gr1_gr2 <- union(gr1$var_id, gr2$var_id)
  
  # Vector of sample-independent annotations
  #-----------------------------------------
  
  stable_annotations <- c("var_id",
                          "Gene.refGene","Func.refGene","ExonicFunc.refGene",
                          "CLNSIG","CLNREVSTAT","CLNDN",
                          "avsnp150","HGVS_cDNA","HGVS_Prot",
                          "Ref","Alt",
                          "CDS_Pos","TranscriptSize",
                          "AA_Pos","AA1","AA2","ProteinSize","Interpro_domain",
                          "gnomAD_genome_ALL","ExAC_ALL")
  
  # Matrices for sample-dependent annotations
  #------------------------------------------
  
  # Matrix for genotypes
  gt.mx <- matrix(NA, ncol=2, nrow=length(union_gr1_gr2))
  rownames(gt.mx) <- union_gr1_gr2
  colnames(gt.mx) <- samples.df$id
  
  # Matrix for BASS scores
  score.mx <- matrix(NA, ncol=2, nrow=length(union_gr1_gr2))
  rownames(score.mx) <- union_gr1_gr2
  colnames(score.mx) <- samples.df$id
  
  # Process non-overlapping variants
  #---------------------------------
  
  # GRanges
  unique_gr1.gr <- gr1[unique_gr1, stable_annotations]
  unique_gr2.gr <- gr2[unique_gr2, stable_annotations]
  
  # Matrices
  gt.mx[unique_gr1, samples.df$id[1] ] <- as.data.frame(gr1[unique_gr1,])$GT
  score.mx[unique_gr1, samples.df$id[1] ] <- as.data.frame(gr1[unique_gr1,])$Score
  
  score.mx[unique_gr2, samples.df$id[2] ] <- as.data.frame(gr2[unique_gr2,])$Score
  gt.mx[unique_gr2, samples.df$id[2] ] <- as.data.frame(gr2[unique_gr2,])$GT
  
  # Process overlapping variants
  #-----------------------------
  
  # Check that the stable annotations are the same before making GRange
  intersect_gr1.gr <- gr1[intersect_gr1_gr2, stable_annotations]
  intersect_gr2.gr <- gr2[intersect_gr1_gr2, stable_annotations]
  if (!all(intersect_gr1.gr==intersect_gr2.gr)) {
    stop("GRanges shared elements differ", ok)
  }
  
  # Grange
  intersect_gr1_gr2.gr <- intersect_gr1.gr
  
  # Matrices
  score.mx[intersect_gr1_gr2, samples.df$id[1] ] <- as.data.frame(gr1[intersect_gr1_gr2,])$Score
  gt.mx[intersect_gr1_gr2, samples.df$id[1] ] <- as.data.frame(gr1[intersect_gr1_gr2,])$GT
  
  score.mx[intersect_gr1_gr2, samples.df$id[2] ] <- as.data.frame(gr2[intersect_gr1_gr2,])$Score
  gt.mx[intersect_gr1_gr2, samples.df$id[2] ] <- as.data.frame(gr2[intersect_gr1_gr2,])$GT
  
  # Combine overlapping and non-overlapping granges
  #------------------------------------------------
  
  # Combine granges
  gr1_gr2 <- c(unique_gr1.gr, unique_gr2.gr, intersect_gr1_gr2.gr)
  
  # Sort granges
  gr1_gr2 <- sort(gr1_gr2)
  
  # Update matrices to keep order consistent with GRange
  score.mx <- score.mx[gr1_gr2$var_id,]
  gt.mx <- gt.mx[gr1_gr2$var_id,]
  
  # Make SummarizedExperiment
  #--------------------------
  
  # Sanity checks: make sure that the objects are in sync
  if (any(colnames(gt.mx) != colnames(score.mx))) {
    stop("Matrices have different samples", ok)
  }
  if (any(rownames(gt.mx) != rownames(score.mx))) {
    stop("Matrices have different variants", ok)
  }
  if (any(names(gr1_gr2) != rownames(gt.mx))) {
    stop("Matrices and GRanges have different variants", ok)
  }
  
  # Make some metadata (optional)
  metadata.ls <- list(source="Made from two GRanges")
  
  # Put all data to a SummarizedExperiment container
  merged.se <- SummarizedExperiment(
    assays=list(genotypes=gt.mx, scores=score.mx),
    rowRanges=gr1_gr2, 
    colData=samples.df, 
    metadata=metadata.ls)
  
  # Return result
  #--------------
  
  return(merged.se)
}

# Pregress report
cat("Loaded function combine_two_GRanges.rms()","\n")
