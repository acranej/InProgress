require(data.table)
require(optparse)
require(parallel)

#' Germline filtering
#' @name fuzzy_filter_germline
#' 
#' @param i: passed from lapply to iterate
#' @param bedpe: Manta Bedpe returned from annotate_sv function
#' @return SV data table with columns added indicating germline or somatic, germline is defined as <=200bp away from agnostic perfect match in reference
#' @description Determines if each SV should be considered germline by hard filtering
#' @export
fuzzy_filter_germline <- function(i, bed) {
  sub <- bed[i,]
  ## reorder for filtering
  if(sub$CHROM_A > sub$CHROM_B) {
    sub_ord <- cbind(CHROM_A=sub$CHROM_B, START_A=sub$START_B, END_A=sub$END_B, CHROM_B=sub$CHROM_A, START_B=sub$START_A, END_B=sub$END_A, sub[,7:ncol(bed)])
  } else if (sub$CHROM_A == sub$CHROM_B & sub$START_A > sub$START_B) {
    sub_ord <- cbind(CHROM_A=sub$CHROM_B, START_A=sub$START_B, END_A=sub$END_B, CHROM_B=sub$CHROM_A, START_B=sub$START_A, END_B=sub$END_A, sub[,7:ncol(bed)])
  } else {
    sub_ord <- sub
  }
  ### change to integers to match reference germline
  sub_ord[CHROM_A == "X", CHROM_A := 23]
  sub_ord[CHROM_B == "X", CHROM_B := 23]
  sub_ord[CHROM_A == "Y", CHROM_A := 24]
  sub_ord[CHROM_B == "Y", CHROM_B := 24]
  ### subset reference to matching chromosome
  ref_sub <- hg38_germline_gnomad[chrom1 == sub_ord$CHROM_A & chrom2 == sub_ord$CHROM_B]
  ### calculate distances
  ref_sub[,str_dist := abs(start - sub_ord$START_A)]
  ref_sub[,end_dist := abs(end - sub_ord$START_B)]
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

option_list <- list(
  make_option(c("-i","--input"), type = "character", default = NULL, help = "Input file"),
  make_option(c("-o","--output"), type = "character", default = NULL, help ="Output directory"),
  make_option(c("-g","--germref"), type = "character", default = NULL, help = "Germline reference"),
  make_option(c("-c","--cores"), type = "integer", default = 3, help = "Number of cores")
)
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)
hg38_germline_gnomad <<- fread(paste0(opt$germref))
input_bed = fread(paste0(opt$input))
output_dir <- paste0(opt$output)

bedpe_fuzzy_filtered <- rbindlist(mclapply(1:nrow(input_bed), fuzzy_filter_germline, input_bed, mc.cores = opt$cores))
output_name <- paste0(output_dir, unlist(strsplit(opt$input, "/"))[length(unlist(strsplit(opt$input, "/")))], "_sv_hg38fuzzyfilter.bedpe")
write.table(bedpe_fuzzy_filtered, output_name, sep = '\t', row.names = F, col.names = T, quote = F)