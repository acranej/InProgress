chrom1=chrom2=str_dist=end_dist=gnomad_dist=start=end=NULL
#' gNOMAD dataset to fuzzy filter
#'
#' gNOMAd v2.1 control sites lifted over to hg38
#' @name gnomad_germline_hg38all
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
gnomad_germline_hg38all = fread(system.file('extdata', 'gnomad_germline_hg38all.txt', package = 'InProgress'))
globalVariables("gnomad_germline_hg38all")

#' @name fuzzy_filter_germline
#' @title Distance to closest germline annotator
#' @param itter passed from mclapply to iterate
#' @param bed Bedpe returned from annotate_sv function
#' @return SV data table with columns added indicating germline or somatic, germline is defined as <=1kbp away from agnostic perfect match in reference
#' @description 
#' 
#' Determines if each SV should be considered germline by hard filtering. Used with wrapper function.
#' 
#' @import data.table 
#' @keywords internal
fuzzy_filter_germline = function(itter = NULL, bed = NULL) {
  sub <- bed[itter,]
  ## reorder for filtering
  if(sub$chrom1 > sub$chrom2) {
    sub_ord <- cbind(chrom1=sub$chrom2, start1=sub$start2, end1=sub$end2, chrom2=sub$chrom1, start2=sub$start1, end2=sub$end1, sub[,7:ncol(bed)])
  } else if (sub$chrom1 == sub$chrom2 & sub$start1 > sub$start2) {
    sub_ord <- cbind(chrom1=sub$chrom2, start1=sub$start2, end1=sub$end2, chrom2=sub$chrom1, start2=sub$start1, end2=sub$end1, sub[,7:ncol(bed)])
  } else {
    sub_ord <- sub
  }
  ### change to integers to match reference germline
  sub_ord[chrom1 == "X", chrom1 := 23]
  sub_ord[chrom2 == "X", chrom2 := 23]
  sub_ord[chrom1 == "Y", chrom1 := 24]
  sub_ord[chrom2 == "Y", chrom2 := 24]
  ### subset reference to matching chromosome
  ref_sub <- gnomad_germline_hg38all[chrom1 == sub_ord$chrom1 & chrom2 == sub_ord$chrom2]
  ### calculate distances
  ref_sub[,str_dist := abs(start - as.numeric(as.character(sub_ord$start1)))]
  ref_sub[,end_dist := abs(end - as.numeric(as.character(sub_ord$start2)))]
  ref_sub[,gnomad_dist := (str_dist + end_dist)]
  ### choose closest match
  ref_min <- ref_sub[which.min(ref_sub$gnomad_dist)]
  if(nrow(ref_min) == 0) {
    sub <- cbind(sub, gnomad_dist = "")
    return(sub)
  } else {
    sub <- cbind(sub, gnomad_dist = ref_min$gnomad_dist)
    return(sub)
  }
}

#' @name closest_germline
#' @title Determines distance to nearest germline event
#' @param bp \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} from \link[InProgress]{svaba_vcf2bedpe} or \link[InProgress]{manta_vcf2bedpe}
#' @param cores Number of cores to run on, default is 1
#' @return \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} with a column added for distance to nearest germline event
#' @description 
#' 
#' Uses \href{https://gnomad.broadinstitute.org/downloads#v2-structural-variants}{gnomAD} to annotate the nearest germline event to each structural variant.
#' For more information read \href{https://www.nature.com/articles/s41586-020-2287-8}{gnomAD blog}. Reference is in hg38.
#' 
#' @import data.table
#' @importFrom parallel mclapply 
#' @export
closest_germline = function(bp = NULL, cores = 1) {
  if(is.null(bp)) {
    stop('NULL input')
  }
  cat("Comparing against known germline...")
  annotated_bedpe <- rbindlist(mclapply(1:nrow(bp), fuzzy_filter_germline, bp, mc.cores = cores))
  cat("done.\n")
  return(annotated_bedpe)
}
