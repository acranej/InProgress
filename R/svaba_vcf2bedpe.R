#' SVaBa vcf to bedpe
#' @name svaba_vcf2bedpe
#' @title Converts SVaBa VCFs to bedpe format
#' @param file_path: file path to SVaBa VCF 
#' @return A bedpe filtered from the VCF
#' @import data.table
#' @importFrom parallel mclapply
#' @description Converts SVaBa VCF to bedpe file format
#' @export

svaba_vcf2bedpe = function(filepath = NULL) {
  if (!file.exists(filepath)) {
    print(paste("File does not exist",filepath))
  }
  vcf_dt <- fread(cmd=paste("grep -v '^#'", filepath),sep='\t')
  
  # Set colnames of vcf_dt to standard...
  if (nrow(vcf_dt) == 0) {
    return (vcf_dt)
  }
  if (ncol(vcf_dt)==10) {
    setnames(vcf_dt, paste0("V",seq(1:10)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL"))
  } else {
    setnames(vcf_dt, paste0("V",seq(1:11)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL","TUMOR"), skip_absent=TRUE)
  }
  
  # extract info column
  if ("INFO" %in% colnames(vcf_dt) ) {
    # Extract info from columns... 
    vcf_dt[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    vcf_dt$sample = gsub("(.*?)_.*","\\1",basename(filepath))
    vcf_dt[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    vcf_dt[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, NUMPARTS := as.integer(gsub(".*?NUMPARTS=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, SCTG := gsub(".*?SCTG=(.*?);.*", "\\1", INFO)]
    vcf_dt[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, REPSEQ := gsub(".*?;REPSEQ=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, REPSEQ := ifelse(grepl(";", REPSEQ), "", REPSEQ)] 
    vcf_dt[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    vcf_dt[, HOMLEN := nchar(HOMSEQ)] 
    vcf_dt[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    vcf_dt[, NDISC := as.numeric(gsub(".*?NDISC=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, SVMETHOD := substr(INFO,regexpr("SVMETHOD=",INFO)+nchar("SVMETHOD="),regexpr(";NDISC",INFO)-1)]
    
  }
  ### only take those that pass the filter
  vcf_dt <- vcf_dt[FILTER == 'PASS']
  
  # More extraction regexpr stuff...
  if ("TUMOR" %in% colnames(vcf_dt)) {
    vcf_dt[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid]
    vcf_dt[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid]
    vcf_dt[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid]
  }
  if ("NORMAL" %in% colnames(vcf_dt)) {
    vcf_dt[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid]
    vcf_dt[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid]
    vcf_dt[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid]
  }
  
  vcf_dt[, strand1 := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
  vcf_dt[, strand2 := rev(strand1), by=uid]
  vcf_dt[, start2 := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  vcf_dt[, chrom2 := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
  vcf_dt[, end1 := start]
  
  bad.ix <- vcf_dt[grepl("^G|^M", seqnames), uid]
  vcf_dt <- vcf_dt[!uid %in% bad.ix]
  
  # Set SID == name of file path without the extension...
  vcf_dt[, sid := basename(tools::file_path_sans_ext(filepath))]
  
  vcf_dt[, mates_idx := unlist(strsplit(ID, ":"))[1], by = "ID"]
  vcf_dt[, which_mate := unlist(strsplit(ID, ":"))[2], by = "ID"]
  
  temp_bedpe <- NULL
  removed_bnd <- NULL
  for(i in 1:length(unique(vcf_dt$mates_idx))){
    
    foo <- vcf_dt[mates_idx == unique(vcf_dt$mates_idx)[i]]
    
    if(!(nrow(foo)== 2)) {
      mes <- paste0("Breakpoint ",  unique(vcf_dt$mates_idx)[i], " has incorrect number of mates for ", foo$sample, " It has been removed.")
      # excludes this breakpoint from 
      continue = FALSE
      removed_bnd <- rbind(removed_bnd, foo)
      warning(mes[1])
    } else {
      continue = TRUE
    }
    if(continue) {
      #### build bedpe
      foo1 <- foo[which_mate == 1]
      foo2 <- foo[which_mate == 2]
      
      bedpe_base <- as.data.frame(cbind(foo1$seqnames, foo1$start, foo1$end,
                                        foo2$seqnames, foo2$start, foo2$end))
      colnames(bedpe_base) <- c("chrom1", "start1", "end1", "chrom2","start2","end2")
      
      if(!(foo1$sid == foo2$sid)){
        stop("Multiple samples are being processed, one at a time please...")
      }
      
      bedpe_base <- cbind(bedpe_base, paste0(foo$sample[1],"_", foo$mates_idx[1]))
      colnames(bedpe_base)[7] <- "name"
      
      bedpe_base <- cbind(bedpe_base, foo$QUAL[1])
      colnames(bedpe_base)[8] <- "score"
      
      bedpe_base <- cbind(bedpe_base, foo1$strand1[1])
      bedpe_base <- cbind(bedpe_base, foo2$strand1[1])
      colnames(bedpe_base)[9:10] <- c("strand1", "strand2")
      
      refs_alts <- as.data.frame(cbind(foo1$REF[1], foo1$ALT[1],
                                       foo2$REF[1], foo2$ALT[1]))
      colnames(refs_alts) <- c("REF_1","ALT_1","REF_2","ALT_2")
      bedpe_base <- cbind(bedpe_base, refs_alts)
      bedpe_base <- cbind(bedpe_base, foo1[,c("SPAN", "HOMSEQ", "HOMLEN","INSERTION","NDISC","FILTER","sample", "TUMALT")])
      mapqs <- as.data.frame(cbind(foo1$MAPQ[1], foo2$MAPQ[1]))
      colnames(mapqs) <- c("MAPQ_1","MAPQ_2")
      bedpe_base <- cbind(bedpe_base, mapqs)
      
      
      temp_bedpe <- rbind(temp_bedpe, bedpe_base)
    }
  }
  
  bedpe <- as.data.table(temp_bedpe)
  ### remove those less than 50bp in SPAN
  bedpe_final <- bedpe[as.numeric(as.character(SPAN)) >=50 | as.numeric(as.character(SPAN)) == -1]
  return(bedpe_final)
}

