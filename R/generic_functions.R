#' The function splits GRangesList files of samples into list of GRangesList
#' for individual chromosomes in samples
#'
#' @import GenomeInfoDb
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' @importFrom GenomicFeatures getChromInfoFromUCSC
#' @import data.table
#' @import ggplot2





#' 
#'
#'
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
                           
                           ## Checks if only standard chromosomes are included
                           if ( ! all(seqlevels(svars) %in% standardChromosomes(svars)) ) {
                               svars <- keepStandardChromosomes(svars, pruning.mode = "tidy")
                               svars <- keepSeqlevels(svars,
                                                      value = seqlevelsInUse(svars))
                           }
                           
                           ## splits into chromosomes
                           return(split(svars, seqnames(svars) ) )
                       }, 
                       mc.cores = n_cores),
                       snames)
                       

    

    return(out)
}



#' This function counts number of variants across samples and chromosomes.
#' The function outputs a matrix with samples as rows and chromosomes as
#' columns. Counts can be normalized for chromosome lengths with
#' chr.norm = TRUE argument.
#'
#'
#' 
#' @export

count_chr_variants  <- function(chr_split_grangeslist, chr.norm = TRUE) {

    
    counts_list <- lapply(chr_split_grangeslist,
                          function(sample_chr) {
                              sapply(sample_chr, length)
                          } )
    counts_matrix  <- do.call(rbind, counts_list)

    if (chr.norm) {
        chroms  <- colnames(counts_matrix)

        genome_assembly  <- GenomeInfoDb::genome(chr_split_grangeslist [[1]]  [[1]] ) [[1]]

        chrom_sizes  <-  getChromInfoFromUCSC(genome_assembly)
        chrom_sizes  <-  data.table(chrom_sizes)
        setkey(chrom_sizes, chrom)
        
        chrom_lengths  <-  setNames(chrom_sizes[chroms, length], chroms)

        length_factors  <- chrom_lengths / min(chrom_lengths)

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


#' Plot chromosome-wise counts with KOs as facets, colored conditions
#' and grouped for clones
#'
#' 
#' 
#' @export


plot_chrwise_counts  <- function(chr_wise_counts, KOs, treatment, clones) {

    dt = as.data.frame(chr_wise_counts)

    dt$KO  <-  KOs
    dt$treatment  <-  treatment
    dt$clones  <-  clones

    dt.gg  <-  reshape2::melt(dt, id = c("KO", "treatment", "clones"))
    
    dt.gg$sample  <-  paste(dt.gg$KO, dt.gg$treatment, dt.gg$clones, sep ="_")


    p  <-  ggplot(dt.gg, aes(x = variable, y = value, color = treatment,
                         group = sample) ) +
        geom_line() +
        geom_point() +
        facet_wrap ( ~ KO) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
    return(p)
}


#' Get mutational matrix in trinucleotide context for individual chromsomes
#' 
#'
#' 
#' 
#' @export

mut_matrix_chr = function() {
}
