#' Convert 96 channel mutational profiles to 6 
#'
#' @details
#' Mutational profiles of samples for individual chromosomes are transformed
#' into 6 mutation types. Works on a list of matrices.
#'
#' @param SNV_profiles_list list of matrices, result from
#' \code{\link{mut_matrix_chr}} function.
#' @param perMb if TRUE, variant counts are normalized by chromosome length.
#' default = TRUE
#'
#' @return list of matrices with 6-channel substitution type counts for each
#' sample
#'
#' @examples
#' 
#' SNVs_chr = split_chr(SNVs)
#' SNV_chr_muts = mut_matrix_chr(SNVs_chr, ref_genome)
#' SNV_chr_mut_types = get_mut_types(SNV_chr_muts, perMb = TRUE, "hg38")
#'

#' @export

get_mut_types  <- function( SNV_profiles_list, perMb = TRUE, genome_assembly) {

    sample_mut_types = list()

    sample_names  <-  names(SNV_profiles_list)

    assembly_norm_factors = get_chrom_norm_factors(genome_assembly) 
    for (sample in sample_names) {
        smat  <- SNV_profiles_list[[sample]]
        
        out  <-  vapply( MUT_TYPES, function(x) 
            rowSums(smat[,    substr(colnames(smat), 3, 5)== x,drop=FALSE]),
            numeric(nrow(smat))
            )

        if (perMb) {
            if (missing(genome_assembly) ) {
                stop("perMb=TRUE option has to be passed with genome_assembly value: e.g. hg19, hg38")
            }
            chrom_norms = assembly_norm_factors [rownames(smat)]
            out = out / chrom_norms
        }
        
        sample_mut_types[[sample]] = out
    }
    return(sample_mut_types)
}
