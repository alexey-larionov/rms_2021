# Function to add GRanges to SummarizedExperiment in RMS progect

# Alexey Larionov, 26 Feb 2021

# Summary
#--------

# Takes GRange, SummarizedExperiment and a list
# The list contains the sample descriptions
# Outputs SummarizedExperiment with added GRanges and sample data

# Required libraries
#-------------------

require(GenomicRanges)
require(SummarizedExperiment)

# Function
#---------

add_GR_to_SE.rms <- function(gr, se, sample.ls){
  
  # Make gr key
  #------------

  # Select data for the key
  gr_chr <- as.character(gr@seqnames)
  gr_start <- as.data.frame(gr@ranges)$start
  gr_end <- as.data.frame(gr@ranges)$end
  gr_ref <- gr$Ref
  gr_alt <- gr$Alt
  
  # Make and add key to GRanges
  gr_var_id <- paste(gr_chr, gr_start, gr_end, gr_ref, gr_alt, sep="_")
  gr$var_id <- gr_var_id
  
  # Add key to GRanges names
  names(gr) <- gr_var_id
  
  # Make subsets
  #-------------
  
  unique_gr <- setdiff(gr$var_id, names(se))
  unique_se <- setdiff(names(se), gr$var_id)
  intersect_gr_se <- intersect(gr$var_id, names(se))
  union_gr_se <- union(gr$var_id, names(se))

  # Vector of sample-independent variant annotations
  #-------------------------------------------------
  
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
  
  # Genotypes
  gt.mx <- matrix(NA,
                  ncol=length(se$id)+1,
                  nrow=length(union_gr_se))
  rownames(gt.mx) <- union_gr_se
  colnames(gt.mx) <- c(se$id, sample.ls$id)

  # Score
  score.mx <- matrix(NA,
                     ncol=length(se$id)+1,
                     nrow=length(union_gr_se))
  rownames(score.mx) <- union_gr_se
  colnames(score.mx) <- c(se$id, sample.ls$id)

  # Process non-overlapping variants
  #---------------------------------
  
  # GRanges
  unique_gr.gr <- gr[unique_gr,stable_annotations]
  unique_se.gr <- rowRanges(se)[unique_se,stable_annotations]

  # Matrices

  score.mx[unique_gr, sample.ls$id] <- as.data.frame(gr[unique_gr,])$Score
  gt.mx[unique_gr, sample.ls$id] <- as.data.frame(gr[unique_gr,])$GT
  
  score.mx[unique_se, se$id] <- assay(se,"scores")[unique_se,]
  gt.mx[unique_se, se$id] <- assay(se,"genotypes")[unique_se,]
  
  # Process overlapping variants
  #-----------------------------
  
  # Check that the stable annotations are the same in GR and SE
  intersect_gr.gr <- gr[intersect_gr_se,stable_annotations]
  intersect_se.gr <- rowRanges(se)[intersect_gr_se,stable_annotations]
  if (!all(intersect_gr.gr==intersect_se.gr)) {
    stop("GRanges shared elements differ", ok)
  }

  # Grange
  intersect_gr_se.gr <- intersect_gr.gr
  
  # Matrices
  
  score.mx[intersect_gr_se, sample.ls$id] <- as.data.frame(gr[intersect_gr_se,])$Score
  gt.mx[intersect_gr_se, sample.ls$id] <- as.data.frame(gr[intersect_gr_se,])$GT
  
  score.mx[intersect_gr_se,se$id] <- assay(se,"scores")[intersect_gr_se,]
  gt.mx[intersect_gr_se,se$id] <- assay(se,"genotypes")[intersect_gr_se,]
  
  # Combine overlapping and non-overlapping granges
  #------------------------------------------------
  
  # Combine granges
  gr_se.gr <- c(unique_gr.gr,unique_se.gr,intersect_gr_se.gr)
  
  # Sort granges
  gr_se.gr <- sort(gr_se.gr)
  
  # Update matrices to keep order consistent with GRange
  score.mx <- score.mx[gr_se.gr$var_id,]
  gt.mx <- gt.mx[gr_se.gr$var_id,]
  
  # Make new SummarizedExperiment
  #------------------------------
  
  # Check that objects are in sync
  if (any(colnames(gt.mx) != colnames(score.mx))) {
    stop("Matrices have different samples", ok)
  }
  
  if (any(rownames(gt.mx) != rownames(score.mx))) {
    stop("Matrices have different variants", ok)
  }
  
  if (any(names(gr_se.gr) != rownames(gt.mx))) {
    stop("Matrices and GRanges have different variants", ok)
  }
  
  # Prepare samples info
  samples <- c(se$id,sample.ls$id)
  files <- c(se$file,sample.ls$file)
  study <- c(se$study,sample.ls$study)
  samples.df <- DataFrame(id=samples, file=files ,study=study, row.names=samples)
  samples.df
  
  # Make some metadata (optional)
  metadata.ls <- list(source="Made from GR and SE")
  
  # Put all data to a SummarizedExperiment container
  merged.se <- SummarizedExperiment(assays=list(genotypes=gt.mx,scores=score.mx),
                                    rowRanges=gr_se.gr, colData=samples.df, metadata=metadata.ls)
  
  # Return result
  #--------------
  
  return(merged.se)
  
}

# Pregress report
cat("Loaded function add_GR_to_SE.rms()","\n")
