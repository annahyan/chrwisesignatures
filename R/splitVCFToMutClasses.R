#' A function splitting vcf files into SNVs, multinucleotide substitutions and indels
#'
#' @details Splits variants in GRanges into SNVs 
#' 
#' @import GenomeInfoDb
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#' @import GenomicRanges
#'
#' @examples
#' splitVCFtoMutClasses(vcf_granges)
#'
#' @export


splitVCFToMutClasses <- function(vcffile, n_cores) {
    
    snames <- names(vcffile)
    
    indels <- list()
    SNVs <- list()
    multisubs <- list()
    

    if (missing(n_cores))
    {
      n_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(n_cores)))
            n_cores <- ceiling(detectCores() * 0.75)
        else
          n_cores = 1
    }
    
    for(sname in snames) {
        cat(sname, "\n")

        svars = vcffile[[sname]]

        ## Removing multiallelic sites
        
        ## Checking the number of alternative alleles
        alt_num = elementNROWS(svars$ALT)
        ## Removing multiallelic sites
        svars = svars[alt_num == 1, ]
        ## Getting the length of the alternative 
        alt_len = width(unlist(svars$ALT))

        ## Making the R
        if ( is (svars$ALT, "DNAStringSetList")) {
            svars$ALT  <-  as.character(svars$ALT@unlistData)
        }
        
        
        
        indels[[sname]]  <-  svars[ alt_len != width(svars$REF), ]
        
        ## multinucleotide substitutions 
        
        variant_names  <-  sub( "(.*)_.*", "\\1", names(svars) )
        
        ## exlcude sites with the same starting coordinate
        invalid_sites  <-  variant_names[which(duplicated(variant_names)) ]

        svar_nonindels  <-  svars[ (! variant_names %in% invalid_sites) & 
                                   (width(svars$REF) == 1 & alt_len == 1), ]
        
        chr_split  <-  split(svar_nonindels, seqnames(svar_nonindels) )

        chr_wise_multisubs = mclapply(chr_split, function(chr_ranges) {
                chrom_name = seqlevelsInUse(chr_ranges)
                chr_ranges = GenomicRanges::as.data.frame(chr_ranges)
                
                ## consecutive positions list
                cons_poss_list  <-  split(chr_ranges$start, f = cumsum(c (0, diff(chr_ranges$start) > 1)) )
                multinucl_groups  <-  cons_poss_list [
                    vapply(cons_poss_list,length, numeric(1) )  > 1
                ]

                concat_variants_df  <-  lapply (multinucl_groups, function(x) { 
                    df_rr  <-  chr_ranges[ paste0(chrom_name, ":", x, "_"), ]
                    
                    ret_df = data.frame(
                        seqnames = unique(df_rr$seqnames), 
                        start = min(df_rr$start),
                        end = max(df_rr$end),
                        width = nrow(df_rr),
                        strand = "*",
                        paramRangeID = unique(df_rr$paramRangeID),
                        REF = paste0(df_rr$REF, collapse = ""),
                        ALT = paste0(df_rr$ALT, collapse = ""),
                        QUAL = min(df_rr$QUAL),
                        FILTER = ".")
                    return(ret_df)
                } )

                return(bind_rows(concat_variants_df, .id = NULL))

            }, mc.cores = n_cores)
        
        multisub_units  <-  bind_rows(chr_wise_multisubs, .id = NULL )
        
        multisubs[[sname]]  <-  GRanges(multisub_units[!is.na(multisub_units$start), ])
        SNVs[[sname]]  <-  subsetByOverlaps(svar_nonindels, multisubs[[sname]], invert = TRUE)
        
    }
    multisubs = GRangesList(multisubs)
    SNVs = GRangesList(SNVs)
    
    return(list(indels = indels, multisubs = multisubs, SNVs = SNVs))
}

