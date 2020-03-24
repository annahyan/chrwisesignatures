#' The function splits GRangesList files of samples into list of GRangesList
#' for individual chromosomes in samples
#'
#' @importFrom utils getFromNamespace
#' @import GenomeInfoDb
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' @importFrom GenomicFeatures getChromInfoFromUCSC
#' @import data.table
#' @import ggplot2
#' @import MutationalPatterns
#' @import robCompositions




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
        theme_bw(base_size = 13) +
        theme(axis.text.x = element_text(angle = 90),
              )
    return(p)
}



mut_96_occurrences  <-  getFromNamespace("mut_96_occurrences", "MutationalPatterns")



#' Get mutational matrix in trinucleotide context for individual chromsomes
#' 
#'
#' 
#' 
#' @export




mut_matrix_chr = function( chr_split_grangeslist, ref_genome, num_cores) {
    if (missing(num_cores))
    {
        num_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
            num_cores <- detectCores()
        else
            num_cores = 1
    }

    samples_list <- mclapply (as.list(chr_split_grangeslist), function (sample_chr_vcfs)
    {
        sample_chr_rows   <-  lapply (sample_chr_vcfs, function(chr_vcf) {
            type_context  <-  type_context(chr_vcf, ref_genome)
            row  <- mut_96_occurrences(type_context)
            return(row)
        }
        )

        df = data.frame()

        ## Merge the rows into a dataframe.
        for (row in sample_chr_rows)
        {
            if (class (row) == "try-error") stop (row)
            df = rbind (df, row)
        }
        
        colnames(df) = names(row)
        rownames(df) = names(sample_chr_vcfs)
        
        return(df)
        
    }, mc.cores = num_cores)
    
}




#' Get mutation types from mutational profiles of sample lists
#' 
#'
#' 
#' 
#' @export

get_mut_types  <- function( SNV_profiles ) {

    sample_mut_types = list()

    sample_names  <-  names(SNV_profiles)

    for (sample in sample_names) {
        smat  <- SNV_profiles[[sample]]
        
        out  <-  vapply( MUT_TYPES, function(x) 
            rowSums(smat[,    substr(colnames(smat), 3, 5)== x,drop=FALSE]),
            numeric(nrow(smat))
            )
        sample_mut_types[[sample]] = out
    }
    return(sample_mut_types)
}



#' Plot CoDa PCA
#' 
#'
#' 
#' 
#' @export


plot_coda_pca  <-  function(dt_list, sample_classes) {

    plot_material = do.call( rbind,
                            lapply(names(dt_list),
                                   function(x) {
                                       df = dt_list[[x]];
                                       rownames(df) = paste(x, rownames(df),
                                                            sep = ":");
                                       return(df)
                                   })
                            )

    sample_types = rep(smp_classes, vapply(dt_list, nrow, numeric(1)))


    pca_list = pcaCoDa(plot_material)
    outliers = outCoDa(plot_material)
    print(outliers)
    
    rob_pca_scores = data.frame(pca_list$scores)
    rob_pca_scores$col = sample_types

    rownames(rob_pca_scores) = rownames(plot_material)


    p = ggplot(rob_pca_scores, aes(x = Comp.1, y = Comp.2, color = col)) + 
        geom_point(size = 2) +
        scale_color_brewer(palette = "Set1" ) # 
#        scale_color_manual(values = sample(CEMM_COLORS_ALL, 9 ) ) # +
    ## ggtitle(paste0("Rank:", nmf.rank))


    pca_loadings = as.data.frame(pca_list$loadings)
    rownames(pca_loadings) = colnames(plot_material)

    p = p + geom_segment(data = pca_loadings,
                         aes(x = 0, y = 0, xend = Comp.1 * 2,
                             yend = Comp.2 * 2),
                         size = 1, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "black") + 
        annotate("text", x = pca_loadings$Comp.1 * 2.2,
                 y = pca_loadings$Comp.2 * 2.2, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab("PC1") + ylab("PC2")
    return(p)
}



#' Plot classical PCA
#' 
#'
#' 
#' 
#' @export


plot_pca  <-  function(dt_list, sample_classes) {

    plot_material = do.call( rbind,
                            lapply(names(dt_list),
                                   function(x) {
                                       df = dt_list[[x]];
                                       rownames(df) = paste(x, rownames(df),
                                                            sep = ":");
                                       return(df)
                                   })
                            )

    sample_types = rep(smp_classes, vapply(dt_list, nrow, numeric(1)))


    pca_list = prcomp(plot_material)

    
    pca_scores = data.frame(pca_list$x)
    pca_scores$col = sample_types

    rownames(pca_scores) = rownames(plot_material)


    p = ggplot(pca_scores, aes(x = PC1, y = PC2, color = col)) + 
        geom_point(size = 2) +
        scale_color_brewer(palette = "Set1" ) # 
#        scale_color_manual(values = sample(CEMM_COLORS_ALL, 9 ) ) # +
    ## ggtitle(paste0("Rank:", nmf.rank))


    pca_loadings = as.data.frame(pca_list$rotation)
    rownames(pca_loadings) = colnames(plot_material)

    p = p + geom_segment(data = pca_loadings,
                         aes(x = 0, y = 0, xend = PC1 * 120,
                             yend = PC2 * 120),
                         size = 1, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "black") + 
        annotate("text", x = pca_loadings$PC1 * 122,
                 y = pca_loadings$PC2 * 122, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab ("PC1") + ylab("PC2")


    return(p)

}





