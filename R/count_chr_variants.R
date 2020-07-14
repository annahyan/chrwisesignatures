#' Count variants across samples and chromosomes
#'
#' @details
#' Provides automated variant counting after variants are split by chromosomes.
#' Normalization by chromosome length is possible, however, the function does
#' not accounting for cases when parts of chromosomes are lost.
#'
#' @param chr_split_grangeslist grangeslist object, result from split_chr
#' function
#'
#' @param perMb if TRUE, variant counts are normalized by chromosome length.
#' default = TRUE
#'
#' @return counts matrix
#'
#' @examples
#'
#' SNVs_chr = split_chr(SNVs)
#' SNVs_chr = lapply( SNVs_chr, function(x) {
#'    x[["chrY"]] = NULL
#'    return(x)
#' }
#' )
#' SNV_chr_counts = count_chr_variants(SNVs_chr, perMb = TRUE)
#'
#' @import GenomeInfoDb
#' 
#' @export

count_chr_variants  <- function(chr_split_grangeslist, perMb = TRUE) {

    
    counts_list <- lapply(chr_split_grangeslist,
                          function(sample_chr) {
                              sapply(sample_chr, length)
                          } )
    counts_matrix  <- do.call(rbind, counts_list)

    if (perMb) {
        chroms  <- colnames(counts_matrix)

        genome_assembly <- GenomeInfoDb::genome(chr_split_grangeslist [[1]]  [[1]] ) [[1]]

        length_factors <- get_chrom_norm_factors(genome_assembly)[chroms]

        if ((ncol(counts_matrix) == length(length_factors))  &
            all.equal(colnames(counts_matrix), names(length_factors) ) ) {
                counts_matrix  <- counts_matrix / length_factors
            } else {
                stop( "The chromosome normalization may be wrong. \
Naming may be disrupted.
I'm not gonna do this to your research!
Quitting... ")
            }
    }
    return(counts_matrix)
}
