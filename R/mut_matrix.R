#' MutationalPatterns mut_matrix function with customizable thread counts
#' 
#' @details
#' See \link[MutationalPatterns]{mut_matrix}
#' 
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSGenome reference genome object
#' @param n_cores Number of cores. If not provided, 3/4 of cores will be used.
#' @return 96 mutation count matrix
#' 
#' @export

mut_matrix = function (vcf_list, ref_genome, n_cores) 
{
    df = data.frame()


    if (missing(n_cores)) {
        n_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(n_cores)))
            n_cores <- ceiling(detectCores() * 0.75)
    }

    rows <- mclapply(as.list(vcf_list), function(vcf) {
        type_context = type_context(vcf, ref_genome)
        row = mut_96_occurrences(type_context)
        return(row)
    }, mc.cores = n_cores)
    for (row in rows) {
        if (class(row) == "try-error") 
            stop(row)
        df = rbind(df, row)
    }
    names(df) = names(row)
    row.names(df) = names(vcf_list)
    return(t(df))
}
