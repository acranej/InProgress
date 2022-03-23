
#' Ensembl gene coordinates
#'
#' Gene coordinates hg38
#' @name hg38_ensembl_genelocations
#' @docType data
#' @keywords data internal
#' @format \code{data.table}

hg38_ensembl_genelocations = readRDS(system.file("extdata","hg38_ensembl_genelocations_formatted.rds", package = "InProgress"))

#' Ensembl exon coordinates
#'
#' Exon coordinates hg38
#' @name hg38_ensembl_exonlocations
#' @docType data
#' @keywords data internal
#' @format \code{data.table}

hg38_ensembl_exonlocations = readRDS(system.file("extdata","hg38_ensembl_exonlocations_formatted.rds", package = "InProgress"))

#' Gene overlap
#' @name find_gene_overlap
#' @title SV occurance in genes
#' @param bed_g bedpe row
#' @return \code{GRangesList} of genes
#' @description 
#'
#' Determines if an SV is overlapping any genes
#' 
#' @import data.table
#' @import GenomicRanges
#' @import gUtils
#' @keywords internal
#' 

find_gene_overlap = function(bed_g = NULL) {
  
}

#' Exon overlap
#' @name find_exon_overlap
#' @title SV occurance in exons
#' @param bed_e bedpe row
#' @return \code{GRanges} of exonsss
#' @description 
#'
#' Determines if an SV is overlapping any exons
#' 
#' @import data.table
#' @import GenomicRanges
#' @import gUtils
#' @keywords internal
#' 

find_exon_overlap = function(bed_e = NULL) {
  
}

##' @name impact_annotation
##' @title 
##'
##'

#impact_annotation = function(i, bed) {
  
#}



