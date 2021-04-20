#' Get mutational matrix in trinucleotide context for individual chromsomes
#'
#' @details
#' Equivalent to mut_matrix from MutationalPatterns, works on a list of
#' grangeslist with chr-variants.
#'
#' @param chr_split_grangeslist grangeslist object, result from split_chr
#' function
#' @param ref_genome BSGenome reference genome object
#' @param n_cores Number of cores. If not provided, 3/4 of cores will be used.
#'
#' @return a list of 96 mutation count matrices
#' 
#' @importFrom parallel mclapply
#' 
#' @export

mut_matrix_chr  <-  function( chr_split_grangeslist, ref_genome, n_cores) {
    if (missing(n_cores))
    {
        n_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(n_cores)))
            n_cores <- detectCores()
        else
            n_cores = 1
    }

    samples_list <- mclapply (as.list(chr_split_grangeslist), function (sample_chr_vcfs)
    {
        ## sample_chr_rows   <-  lapply (sample_chr_vcfs, function(chr_vcf) {
        ##     type_context  <-  type_context(chr_vcf, ref_genome)
        ##     row  <- mut_96_occurrences(type_context)
        ##     return(row)
        ## }
        ## )

        ## df = data.frame()

        ## ## Merge the rows into a dataframe.
        ## for (row in sample_chr_rows)
        ## {
        ##     if (class (row) == "try-error") stop (row)
        ##     df = rbind (df, row)
        ## }
        
        ## colnames(df) = names(row)
        ## rownames(df) = names(sample_chr_vcfs)
        
        ## return(df)

        mut_matrix(sample_chr_vcfs, ref_genome)
        
    }, mc.cores = n_cores)
    
}
