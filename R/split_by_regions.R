#' Split variants into genomic regions: exons, introns, utr3, utr5, intergenic
#'
#' @details Variants are split based on genomic regions - exons, introns, 3'UTR,
#' 5'UTR, transcripts, intergenic regions defined from TxDb.
#'
#' @param variants_list List of GRanges objects. Each element of a list is the
#' set of sample variants.
#' @param txdb_file TxDb object matching with the genome assembly of variants
#' @param n_cores Number of cores. If not provided, 3/4 of cores will be used.#' 
#'
#' @return A list of genomic regions (exons, introns, etc.), containing a list of
#' samples with variants residing in corresponding regions.
#'
#' @examples
#' SNVs_hg38_split = split_by_regions(SNVs, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' @importFrom parallel mclapply
#' @importFrom IRanges subsetByOverlaps
#' 
#' @export



split_by_regions <- function(variants_list, txdb_file, n_cores) {


    
    if (missing(n_cores)) {
        n_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(n_cores)))
            n_cores <- ceiling(detectCores() * 0.75)
        else
            n_cores = 1
    }

    sample_ids = names(variants_list)

    gr_exons = exons(txdb_file)
    gr_introns = unique(unlist(intronsByTranscript(txdb_file)))
    gr_utr3 = unique(unlist(threeUTRsByTranscript(txdb_file)))
    gr_utr5 = unique(unlist(fiveUTRsByTranscript(txdb_file)))
    gr_transcripts = transcripts(txdb_file)

    regions = c("exons", "introns", "utr3", "utr5", "transcripts", "intergenic")

    variants_by_region = list()

    for (region in regions) {

        variants_by_region[[region]] = list()

        if (region == "intergenic") {
            region_fetch = "transcripts"   
        } else {
            region_fetch = region
        }
        
        region_granges = get(paste0("gr_", region_fetch))

        variants_by_region[[region]] = mclapply(
            variants_list, function(sample_variants) {
                if (region == "intergenic") {
                    out = subsetByOverlaps(sample_variants, region_granges,
                                           invert = TRUE)
                } else {
                    out = subsetByOverlaps(sample_variants, region_granges,
                                           invert = FALSE)
                }
                ## splits into chromosomes
                return(out)
            }, mc.cores = n_cores)
    }

    return(variants_by_region)
}
