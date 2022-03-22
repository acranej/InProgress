str1_dist_l=start1_dist_s=str2_dist_l=start2_dist_s=line_dist=sine_dist=NULL
#' Line repeatmaskers
#'
#' repeat masker line elements in hg38
#' @name LINE_dt
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
LINE_dt = readRDS(system.file('extdata', 'repeatmasker_hg38_LINE.bed', package = 'InProgress'))


#' Sine repeatmaskers
#'
#' repeat masker sine elements in hg38
#' @name SINE_dt
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
SINE_dt = readRDS(system.file('extdata', 'repeatmasker_hg38_SINE.bed', package = 'InProgress'))

#' Check format of bedpe
#' @name check_reformat
#' @title reformat for matching
#' @param i iteration value
#' @param df bedpe to be passed 
#' @return data.table that is ordered with lower position first
#' @description Reorders foor acccurate distance calculation
#' @import data.table
#' @keywords internal
#' 
check_reformat = function(i, df){
  
  row <- df[i,]
  
  if (row$chrom1 > row$chrom2) {
      return(cbind(chrom1=row$chrom2, start1=row$start2, end1=row$start2, chrom2=row$chrom1, 
                 start2=row$start1, end2=row$start1, name=row$name, score=row$score, 
                 strand1=row$strand2, strand2=row$strand1, row[,11:ncol(df)]))
    } else if (row$chrom1==row$chrom2 & row$start1>row$start2){
      return(cbind(chrom1=row$chrom2, start1=row$start2, end1=row$start2, chrom2=row$chrom1, 
                 start2=row$start1, end2=row$start1, name=row$name, score=row$score, 
                 strand1=row$strand2, strand2=row$strand1, row[, 11:ncol(df)]))
    } else {
      return(row)
    }
}

#' Closest line element
#' @name find_closest_match_line
#' @title annotates the closest line element
#' @param i iteration value
#' @param bedpe_l ordered bedpe from \link[InProgress]{check_reformat}
#' @return data.table with closest line element distance
#' @description Annotates the distance to closest line element
#' @import data.table
#' @keywords internal
#' 
find_closest_match_line = function(i, bedpe_l){
  row_l <- bedpe_l[i,]
  
  ### both bedpe and ref should be sorted so lower bkpt comes first 
  ref_sub <- LINE_dt[ seqnames == row_l$chrom1 & seqnames == row_l$chrom2]
  ref_sub[,str1_dist_l := abs(start - as.numeric(as.character(row_l$start1)))]
  ref_sub[,str2_dist_l := abs(end - as.numeric(as.character(row_l$start2)))]
  ref_sub[,line_dist := (str1_dist_l + str2_dist_l)]
  
  
  line_min<-ref_sub[which.min(ref_sub$line_dist),]
  
  if(nrow(line_min) == 0) {
    row_l$line_dist <- ""
    return(row_l)
  } else {
    return(cbind(row_l, line_min$line_dist))
  }
}

#' Closest sine element
#' @name find_closest_match_sine
#' @title annotates the closest sine element
#' @param i iteration value
#' @param bedpe_s ordered bedpe from \link[InProgress]{check_reformat}
#' @return data.table with closest sine element distance
#' @description Annotates the distance to closest sine element
#' @import data.table
#' @keywords internal
#' 
find_closest_match_sine = function(i, bedpe_s){
  row_s <- bedpe_s[i,]
  
  ### both bedpe and ref should be sorted so lower bkpt comes first 
  ref_sub <- SINE_dt[ seqnames == row_s$chrom1 & seqnames == row_s$chrom2]
  ref_sub[,str1_dist_s := abs(start - as.numeric(as.character(row_s$start1)))]
  ref_sub[,str2_dist_s:= abs(end - as.numeric(as.character(row_s$start2)))]
  ref_sub[,sine_dist := (str1_dist_s + str2_dist_s)]
  
  
  sine_min<-ref_sub[which.min(ref_sub$sine_dist),]
  
  if(nrow(sine_min) == 0) {
    row_s$sine_dist <- ""
    return(row_s)
  } else {
    return(cbind(row_s, sine_min$sine_dist))
  }
}

#' Wrapper for LINE/SINE annotation
#' @name closest_line_sine 
#' @title Annotate LINE and SINE elements
#' @param bp \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} from \link[InProgress]{svaba_vcf2bedpe} or \link[InProgress]{manta_vcf2bedpe}
#' @param cores Number of cores to run on, default is 1
#' @return \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} with a columns added for distance to nearest LINE and SINE element
#' @description
#' 
#' Annotates the distance to nearest LINE and SINE element using repeat masker elements from \href{http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html}{RepeatMasker}.
#' Annotation is in hg38.
#' 
#' @import data.table
#' @importFrom parallel mclapply
#' @export
closest_line_sine = function(bp = NULL, cores = 1) {
  
  if(is.null(bp)) {
    stop('NULL input')
  }
  
  cat("Checking format...")
  bp_ord <- rbindlist(mclapply(1:nrow(bp), check_reformat, bp, mc.cores = cores))
  cat("done.\n")
  
  cat("Comparing against LINE elements...")
  line_annotated <- rbindlist(mclapply(1:nrow(bp_ord), find_closest_match_line, bp_ord, mc.cores = cores), use.names = FALSE)
  cat("done.\n")
  
  cat("Comparing against SINE elements...")
  sine_line_annotated <- rbindlist(mclapply(1:nrow(line_annotated), find_closest_match_sine, line_annotated, mc.cores = cores), use.names = FALSE)
  cat("done.\n")
  
  return(sine_line_annotated)
}
