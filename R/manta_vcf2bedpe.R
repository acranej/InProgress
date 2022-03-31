MATE_ID=start1=end2=svtype=chrom1=chrom2=start2=end1=strand1=strand2=REF_1=REF_2=REF=name=NAME_2=ALT_1=ALT_2=INFO_1=INFO_2=uuid=NULL
#' Read mate pairing for bedpe construction, internally used
#' @name bnd_matching
#' @title match breakends
#' @param id mate id to pair breakpoints
#' @param bnd data table of breakend SVs
#' @return SV data table with matched breakpoints
#' @description Pairs mate ids from vcf and constructs bedpe row for the SV
#' @import data.table
#' @keywords internal
#'
bnd_matching = function(id, bnd) {
  which_mate_A <- bnd[grepl(id, name)]
  which_mate_B <- bnd[which_mate_A$name == MATE_ID]
  
  # checks for exactly one pair of breakends to match
  # some SVs have one breakpoint that does not pass filter so those SVs will not result in the final bedpe
  if(nrow(which_mate_A) == 1 && nrow(which_mate_B) == 1) {
    bed_temp <- as.data.table(cbind(which_mate_A$seqnames, which_mate_A$start, which_mate_A$start, 
                                    which_mate_B$seqnames, which_mate_B$start, which_mate_B$start))
    colnames(bed_temp) <- c("chrom1","start1","end1","chrom2","start2","end2")
    bed_temp[is.na(end1), end1 := start1]
    bed_temp[is.na(end2), end2 := start2]
    
    #### determine orientation
    strA <- strsplit(which_mate_A$ALT, "")
    str_m <- grep("[[]",strA)
    str_p <- grep("[]]", strA)
    
    ### determines +/- of strand B
    strandB <- ifelse((length(str_m) > 0), "-", ifelse((length(str_p) > 0), "+", "*"))
    
    ### determine orientation
    strB <- strsplit(which_mate_B$ALT, "")
    str_m <- grep("[[]",strB)
    str_p <- grep("[]]", strB)
    
    ### determines +/- of strand A
    strandA <- ifelse((length(str_m) > 0), "-", ifelse((length(str_p) > 0), "+", "*"))
    
    bed_temp <- cbind(bed_temp, which_mate_A$name, which_mate_A$score, strandA, strandB, which_mate_A$svtype,
                      which_mate_A$FILTER, which_mate_A$REF, 
                      which_mate_A$ALT, which_mate_A$MATE_ID, which_mate_B$REF, which_mate_B$ALT,
                      which_mate_A$INFO, which_mate_B$INFO, which_mate_A$FORMAT, which_mate_A$OCILY12, which_mate_A$SPAN)
    
    colnames(bed_temp)[7:22] <- c("name", "score", "strand1","strand2","svtype",
                                  "FILTER","REF_1", 
                                  "ALT_1","MATE_ID","REF_2","ALT_2",
                                  "INFO_1","INFO_2", "FORMAT", "OCILY12","SPAN")
    
    bed_temp <- cbind(bed_temp, which_mate_A$HOMSEQ, which_mate_A$HOMLEN)
    colnames(bed_temp)[23:24] <- c("HOMSEQ", "HOMLEN")
    bed_temp[, uuid := paste0(name, ":",chrom1,
                              ":", start1,":", end1,":", chrom2, ":", start2,
                              ":", end2)]
    
    
    
    return(bed_temp)
  }
}

#' 
#' Function to call for vcf to bedpe conversion
#' @name manta_vcf2bedpe
#' @title Converts \href{https://github.com/Illumina/manta}{MANTA} VCFs to \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} format
#' @param filepath file path to \href{https://github.com/Illumina/manta}{MANTA} vcf
#' @return SV data table in \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} format
#' @description Creates a filtered \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} file from an input \href{https://github.com/Illumina/manta}{MANTA} vcf
#' @import data.table
#' @export
#'
manta_vcf2bedpe = function(filepath= NULL) {
  cat("Reading manta vcf...")
  
  vcf.input <- data.table::fread(cmd=paste("grep -v '^#'", filepath), sep='\t')
  
  if (nrow(vcf.input) == 0) {stop(sprintf('This file is empty!'))}
  
  if (ncol(vcf.input)>=10) {
    setnames(vcf.input, paste0("V",seq(1:10)), c("seqnames","start","name","REF",
                                                 "ALT","score","FILTER","INFO",
                                                 "FORMAT","OCILY12"))
  }
  
  if ("INFO" %in% colnames(vcf.input)) {
    # Extract info from columns...
    
    ### SV TYPE
    vcf.input[, svtype := gsub(".*?SVTYPE=([A-Z]+).*", "\\1", INFO)]
    
    ### SPAN
    vcf.input[, SPAN := as.numeric(gsub(".*?SVLEN=([-0-9]+).*","\\1",INFO))]
    vcf.input$SPAN <-  gsub("-","", vcf.input$SPAN)
    vcf.input[svtype == "BND", SPAN := -1] 
    vcf.input$seqnames <- gsub("chr", "",vcf.input$seqnames)
    
    ### remove those that are not over 50 bp
    
    vcf.input_s <- vcf.input[as.numeric(SPAN) >= 50 | as.numeric(SPAN) == -1]
    
    ### remove imprecise calls 
    
    vcf.input_ss <- vcf.input_s[!(grepl("IMPRECISE", INFO))]
    
    ### HOMSEQ
    vcf.input_ss[grepl("HOMSEQ", INFO), HOMSEQ := gsub(".*?HOMSEQ=([A-Z]+).*","\\1",INFO)]
    
    ### HOMLEN
    vcf.input_ss[, HOMLEN := as.numeric(gsub(".*?HOMLEN=([-0-9]+).*","\\1",INFO))]
    vcf.input_ss[is.na(HOMLEN), HOMSEQ := ""]
    
    ### remove non standard seqs (chromosomes), also for alts
    vcf.input_com_chrom <- vcf.input_ss[seqnames %in% c(1:22,"X", "Y")] 
    vcf.input_com_chrom_s <- vcf.input_com_chrom[!(grepl("chrUn",ALT))]
    vcf.input_com_chrom_ss <- vcf.input_com_chrom_s[!(grepl("random",ALT))]
    vcf.input_com_chrom_sss <- vcf.input_com_chrom_ss[!(grepl("alt",ALT))]

    
    #### build bedpe from non BND calls
    non_bnd <- vcf.input_com_chrom_sss[svtype %in% c("DEL","DUP","INV")]
    
    ### renaming chrom 
    non_bnd[,chrom1 := seqnames ]
    non_bnd[,chrom2 := seqnames ]
    
    ### creating point A positions
    non_bnd[, start1 := start]
    non_bnd[, end1 := start1]
    non_bnd[is.na(end1), end1 := start1]
    
    ### creating point B positions
    non_bnd[, start2 := as.numeric(gsub(".*?END=([-0-9]+).*","\\1",INFO))]
    non_bnd[, end2 := start2]
    non_bnd[is.na(end2), end2 := start2]
    
    ### strand orientations  ---> need to figure out a better way to do this, theese are fake
    non_bnd[, strand1 := "+"]
    non_bnd[, strand2 := "-"]
    
    ### ref and alt poiint A
    non_bnd[, REF_1 := REF]
    non_bnd[, ALT_1 := ALT]
    
    ### ref and alt poiint B
    non_bnd[, REF_2 := "."]
    non_bnd[, ALT_2 := "."]
    
    non_bnd[, INFO_1 := INFO]
    non_bnd[, INFO_2 := "."]
    
    ### unique ID for each SV in thee file
    non_bnd[, uuid := paste0(name, ":",chrom1,
                             ":", start1,":", end1,":", chrom2, ":", start2,
                             ":", end2)]
    
    ### subseting to processeed data table
    non_bnd_save <- non_bnd[,c("chrom1","start1","end1","chrom2","start2",
                               "end2", "name", "score", "strand1","strand2","svtype",
                               "FILTER","REF_1", "ALT_1","REF_2",
                               "ALT_2","INFO_1","INFO_2", "FORMAT", "OCILY12","SPAN",
                               "HOMSEQ", "HOMLEN","uuid")]
    
    
    #### build bedpe from BND
    ### translocations must be formatted separately and must have their mate ids matched since they are located on different chromosomes 
    bnd_ <- vcf.input_com_chrom_sss[svtype == "BND"]
    
    bnd_[,MATE_ID := unlist(strsplit(unlist(strsplit(INFO, "MATEID="))[2],"[;]"))[1], by = "name"]
    cat("done.\n")
    cat('Building bedpe...')
    bnd_bed <- lapply(bnd_$name, bnd_matching, bnd = bnd_)
    bnd_bed <- rbindlist(bnd_bed)
    
    ### prevents duplication of translocation SVs in bedpe with the same positions
    which_dup <- NULL
    for(i in 1:nrow(bnd_bed)) {
      tmp_id <- bnd_bed$name[i]
      id_mate <- which(bnd_bed$MATE_ID == tmp_id)
      if (i < id_mate) {
        which_dup <- rbind(which_dup, cbind(i, id_mate))
      } else {
        which_dup <- rbind(which_dup, cbind(id_mate, i))
      }
    }
    which_dup_dedup <- as.data.frame(which_dup[!(duplicated(which_dup)),])
    
    bnd_bed_dedup <- bnd_bed[which_dup_dedup$i,]
    bnd_bed_dedup[,MATE_ID := NULL]
    ### creates final bedpe and sets positions as numeric
    bedpe <- rbind(non_bnd_save, bnd_bed_dedup)
    bedpe$start1 <- as.numeric(as.character(bedpe$start1))
    bedpe$start2 <- as.numeric(as.character(bedpe$start2))
    bedpe$end1 <- as.numeric(as.character(bedpe$end1))
    bedpe$end2 <- as.numeric(as.character(bedpe$end2))
  } else {
    stop("VCF format looks strange, cannot build bedpe.")
  }
  cat("done.\n")
  return(bedpe)
}
