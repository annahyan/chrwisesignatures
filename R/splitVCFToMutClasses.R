#' A function splitting vcf files into SNVs, multinucleotide substitutions and indels
#'
#' @examples
#' splitVCFtoMutClasses(vcf_granges)
#'
#' @export


splitVCFToMutClasses = function(vcffile) {
  
  snames <- names(vcffile)
  
  indels <- list()
  SNVs <- list()
  multisubs <- list()
  
  for(sname in snames) {
    cat(sname, "\n")

    ## Removing non-canonical chromsomes in sample variants
    svars <- unique(vcffile[[sname]])
    svars <- GenomeInfoDb::keepStandardChromosomes(svars, pruning.mode = "tidy")
    svars <- GenomeInfoDb::keepSeqlevels(svars,
                                         value = GenomicInfoDb::seqlevelsInUse(svars))
    
    ## Removing multiallelic sites
    
    alt_len  <-  sapply(svars$ALT, function(x) { 
        if (length(x) > 1) {
            return(0)
        } else {
            return ( nchar( x[[1]] ) )
        } } 
        )
    
    svars  <-  svars[alt_len > 0, ]
    alt_len  <-  alt_len [alt_len > 0 ]
    
    
    svars$ALT  <-  as.character(svars$ALT@unlistData)
    
    indels[[sname]]  <-  svars[ alt_len != nchar(svars$REF), ]
    
    ## multinucleotide substitutions 
  
    variant_names  <-  sub( "(.*)_.*", "\\1", names(svars) )
        
    ## exlcude sites with the same starting coordinate
    invalid_sites  <-  variant.names[which(duplicated(variant_names)) ]

    svar_nonindels  <-  svars[ (! variant_names %in% invalid_sites) & 
                               (nchar(svars$REF) == 1 & alt_len == 1), ]
    
    chr_split  <-  split(svar_nonindels, seqnames(svar_nonindels) )
    
    multisub_units  <-  do.call( 
        rbind, lapply(chr_split, function(chr_ranges) {
            chrom_name = GenomeInfoDb::seqlevelsInUse(chr_ranges)
            chr_ranges = GenomicRanges::as.data.frame(chr_ranges)
            
            
            ## consecutive positions list
            cons_poss_list  <-  split(chr_ranges$start, f = cumsum(c (0, diff(chr.ranges$start) > 1)) )
            multinucl_groups  <-  cons.poss.list[ vapply(cons_poss_list, length) > 1 ]
            
            do.call(rbind, lapply ( multinucl_groups, function(x) { 
                df_rr  <-  chr_ranges[ paste0(chrom_name, ":", x, "_"), ]
                
                ret_df = data.frame(
                    seqnames = unique(df.rr$seqnames), 
                    start = min(df.rr$start),
                    end = max(df.rr$end),
                    width = nrow(df.rr),
                    strand = "*",
                    paramRangeID = unique(df.rr$paramRangeID),
                    REF = paste0(df.rr$REF, collapse = ""),
                    ALT = paste0(df.rr$ALT, collapse = ""),
                    QUAL = min(df.rr$QUAL),
                    FILTER = ".")
                return(df_rr)
            } ) 
            )
        } ) )
    
    multisubs[[sname]]  <-  GRanges(multisub_units)
    SNVs[[sname]]  <-  subsetByOverlaps(svar_nonindels, multisubs[[sname]], invert = TRUE)
    
  }
    return(list(indels = indels, multisubs = multisubs, SNVs = SNVs))
}


