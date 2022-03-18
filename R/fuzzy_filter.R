chrom1=chrom2=str_dist=end_dist=tot_dist=start=end=NULL
#' gNOMAD dataset to fuzzy filter
#'
#' gNOMAd v2.1 control sites lifted over to hg38
#' @name gnomad_germline_hg38all
#' @docType data
#' @keywords data
#' @format \code{data.table}
hg38_germline_gnomad = fread(system.file('extdata', 'gnomad_germline_hg38all.txt', package = 'InProgress'))
globalVariables("hg38_germline_gnomad")

#' @name fuzzy_filter_germline
#' @title Distance to closest germline annotator
#' @param itter passed from mclapply to iterate
#' @param bed Bedpe returned from annotate_sv function
#' @return SV data table with columns added indicating germline or somatic, germline is defined as <=1kbp away from agnostic perfect match in reference
#' @description 
#' 
#' Determines if each SV should be considered germline by hard filtering, used with wrapper function
#' 
#' @import data.table 
#' @importFrom parallel mclapply
#' @export

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
  ref_sub <- hg38_germline_gnomad[chrom1 == sub_ord$chrom1 & chrom2 == sub_ord$chrom2]
  ### calculate distances
  ref_sub[,str_dist := abs(start - sub_ord$start1)]
  ref_sub[,end_dist := abs(end - sub_ord$start2)]
  ref_sub[, tot_dist := (str_dist + end_dist)]
  ### choose closest match
  ref_min <- ref_sub[which.min(ref_sub$tot_dist)]
  sub <- cbind(sub, ref_min[,c(7:9)])
  if (nrow(ref_min) < 1) {
    sub$Filter[1] <- "Somatic"
  } else if (ref_min$tot_dist > 1000 & sub$SPAN > 0 & sub$SPAN < 1000) {
    sub$Filter[1] <- "Germline"
  } else if (ref_min$tot_dist > 1000) {
    sub$Filter[1] <-paste0("Somatic")
  } else {
    sub$Filter[1] <-paste0("Germline")
  }
  return(sub)
}
