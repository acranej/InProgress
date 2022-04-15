genelocations=exonlocations=NULL
#' Gencode gene coordinates, subset to protein coding
#'
#' Gene coordinates hg38
#' @name hg38_gencode_genelocations
#' @docType data
#' @keywords data internal
#' @format \code{data.table}

hg38_gencode_genelocations = readRDS(system.file("extdata","hg38_geneRanges.rds", package = "InProgress"))

#' Gencode gene coordinates, subset to protein coding
#'
#' Gene coordinates hg19 
#' @name hg19_gencode_genelocations
#' @docType data
#' @keywords data internal
#' @format \code{data.table}

hg19_gencode_genelocations = readRDS(system.file("extdata","hg19_geneRanges.rds", package = "InProgress"))

#' Ensembl exon coordinates
#'
#' Exon coordinates hg38
#' @name hg38_ensembl_exonlocations
#' @docType data
#' @keywords data internal
#' @format \code{data.table}

hg38_ensembl_exonlocations = readRDS(system.file("extdata","hg38_exonRanges.rds", package = "InProgress"))

#' Ensembl exon coordinates
#'
#' Exon coordinates hg19
#' @name hg19_exonlocations
#' @docType data
#' @keywords data internal
#' @format \code{data.table}

hg19_exonlocations = readRDS(system.file("extdata","hg19_exonRanges.rds", package = "InProgress"))

#' @name impact_annotation
#' @title Clqssificatiion decision tree for annotation
#' @param i itteration value
#' @param bed file to be itterated over
#' @return row with annotation column
#' @description 
#' 
#' Determines if an SV is Translocation, Non-protein coding, intron, CN, or coding
#' @import gUtils
#' @import IRanges
#' @import GenomicRanges
#' @keywords internal

impact_annotation = function(i, bed, genelocations = NULL, exonlocations = NULL) {
  sub <- bed[i,]
  
  ### TRA if SV is translocations
  if(sub$chrom1 != sub$chrom2) {
    sub <- cbind(sub, funct_annot = "TRA")
    return(sub)
  }
  
  ### if its intrachromosomal, compare to genes
  sub_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start1)), as.numeric(as.character(sub$start2))))
  bp1_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start1)), width = 1))
  bp2_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start2)), width = 1))
  
  ### first see if it is in a gene
  genes_overlapped <- genelocations %&% sub_gr
  
  bp1_gene <- genelocations %&% bp1_gr
  bp2_gene <- genelocations %&% bp2_gr
  bp1_bp2_gene_count <- length(bp1_gene) + length(bp2_gene)
  
  ### if there is no gene overlap, there will be no exon overlap, so this is non_protein coding
  if(length(genes_overlapped) == 0) {
    sub <- cbind(sub, funct_annot = "Non_Coding")
    return(sub)
  }
  
  ### if it overlaps more than 1  gene and is del/dup, then it has to have a CN impact 
  if(length(genes_overlapped) > 0 & sub$svtype %in% c('DEL','DUP') & bp1_bp2_gene_count == 0) {
    sub <- cbind(sub, funct_annot = "CN")
    return(sub)
  }
  
  ### INV?? --> think about if this actually useful, kinda neutral event except regulatory elements?? idk, will ponder more
  if(length(genes_overlapped) > 0 & sub$svtype %in% c('h2hINV', 't2tINV','INV') & bp1_bp2_gene_count == 0) {
    sub <- cbind(sub, funct_annot = "INV")
    return(sub)
  }
  
  #### those that have bp in genes !!!
  if(bp1_bp2_gene_count > 0) {
    exons_bp1_gene <-exonlocations %&% bp1_gene
    exons_overlapped_bp1 <- exons_bp1_gene %&% sub_gr
    
    exons_bp2_gene <- exonlocations %&% bp2_gene
    exons_overlapped_bp2 <- exons_bp2_gene %&% sub_gr
    
    exons_impacted_in_bp_genes <- length(exons_overlapped_bp1) + length(exons_overlapped_bp2)
    
    # if the breakpoint impacts exons in a gene it has a coding effect
    if(exons_impacted_in_bp_genes > 0) {
      sub <- cbind(sub, funct_annot = "Coding")
      return(sub)
    }
    
    ## if it impacts no exons and is in only 1 gene it is intronic
    if(exons_impacted_in_bp_genes == 0 & length(genes_overlapped) > 0) {
      sub <- cbind(sub, funct_annot = "Intron")
      return(sub)
    } else if(length(genes_overlapped) == 1 & sub$svtype %in% c('DEL','DUP') & bp1_bp2_gene_count == 0) {
      sub <- cbind(sub, funct_annot = "CN")
      return(sub)
    } else if(length(genes_overlapped) == 1 & sub$svtype %in% c('h2hINV', 't2tINV','INV') & bp1_bp2_gene_count == 0) {
      sub <- cbind(sub, funct_annot = "INV")
      return(sub)
    } else {
      sub <- cbind(sub, funct_annot = "")
      return(sub)
    }
  }
}

#' @name funct_annot 
#' @title Annotates SVs
#' @param bp bedpe file to annotate
#' @param cores number of cores to run on, default is 1
#' @param genome run under hg19 or hg38
#' @description 
#' 
#' Annotates for introns, coding, cn, inv, non coding, sv impacts
#'
#' @return bedpe with funct_annot column
#' @import data.table
#' @importFrom parallel mclapply
#' @export


funct_annot = function(bp = NULL, genome = NULL, cores = 1) {
  if(is.null(bp)) {
    stop('NULL input')
  }
  
  if(as.character(genome) == 'hg19') {
    genelocations = hg19_gencode_genelocations
    exonlocations = hg19_exonlocations
  } else if(as.character(genome) == 'hg38') {
    genelocations = hg38_gencode_genelocations
    exonlocations = hg38_ensembl_exonlocations
  } else {
    stop('Please state hg19 or hg38 as genome')
  }
  cat("Annotating SVs...")
  funct_annot_bp <- rbindlist(mclapply(1:nrow(bp), impact_annotation, bp, genelocations = genelocations, exonlocations = exonlocations, mc.cores = cores))
  cat("done.\n")
  return(funct_annot_bp)
}

