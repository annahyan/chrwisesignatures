#' The function splits GRangesList files of samples into into GRangesList
#' for individual chromosomes
#'
#' import GenomeInfoDb
#' import parallel detectCores
#' import parallel mclapply
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
    
    out  <- mclapply(snames, function(sname) {

        svars = vcf_granges_list[[sname]]
        
        ## Checks if only standard chromosomes are included
        if ( ! all(seqlevels(svars) %in% standardChromosomes(svars)) ) {
            svars <- keepStandardChromosomes(svars, pruning.mode = "tidy")
            svars <- keepSeqlevels(svars,
                                   value = seqlevelsInUse(svars))
        }

        ## splits into chromosomes
        return(split(svars, seqlevels(svars)))
    }, 
    mc.cores = n_cores)

    return(GRangesList(out))
}
