#' Split variants into chromosomes
#'
#' @details
#' The function splits GRangesList objects of samples into a list of GRangesList
#' with variants from individual chromosomes.
#'
#' @param vcf_granges_list GRangesList GRanges variants across samples
#' @param n_cores Number of cores. If not provided, 3/4 of cores will be used.
#'
#' @return Returns a list of GRangesLists, for each sample the variants are
#' split by chromosomes.
#'
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' 
#' 
#' @export


split_chr <- function(vcf_granges_list, n_cores) {

    snames <- names(vcf_granges_list)

    if (missing(n_cores)) {
        n_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(n_cores)))
            n_cores <- ceiling(detectCores() * 0.75)
        else
            n_cores = 1
    }
    
    out  <- stats::setNames (
                       
                       mclapply(snames, function(sname) {
                           
                           svars = vcf_granges_list[[sname]]
                           ## splits into chromosomes
                           return(split(svars, seqnames(svars) ) )
                       }, 
                       mc.cores = n_cores),
                       snames)
    return(out)
}
