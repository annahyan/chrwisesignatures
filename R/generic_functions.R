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
#' @importFrom ggpubr ggarrange



get_chrom_norm_factors = function(genome_assembly) {

    chrom_sizes  <-  getChromInfoFromUCSC(genome_assembly)
    chrom_sizes  <-  data.table(chrom_sizes)
    setkey(chrom_sizes, chrom)
    
    ## chrom_lengths  <-  setNames(chrom_sizes[chroms, length], chroms)

    chrom_lengths = setNames(chrom_sizes$length, chrom_sizes$chrom)
    ## length_factors  <- chrom_lengths / min(chrom_lengths)

    length_factors  <- chrom_lengths / 10^6

    return(length_factors)
}

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
#' perMb = TRUE argument.
#'
#'
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

        genome_assembly  <- GenomeInfoDb::genome(chr_split_grangeslist [[1]]  [[1]] ) [[1]]

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
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("Chromosomes") + ylab("Counts")
    
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

get_mut_types  <- function( SNV_profiles, perMb = FALSE, genome_assembly) {

    sample_mut_types = list()

    sample_names  <-  names(SNV_profiles)

    for (sample in sample_names) {
        smat  <- SNV_profiles[[sample]]
        
        out  <-  vapply( MUT_TYPES, function(x) 
            rowSums(smat[,    substr(colnames(smat), 3, 5)== x,drop=FALSE]),
            numeric(nrow(smat))
            )

        if (perMb) {
            if (missing(genome_assembly) ) {
                stop("perMb=TRUE option has to be passed with genome_assembly value: e.g. hg19, hg38")
            }
            chrom_norms = get_chrom_norm_factors(genome_assembly) [rownames(smat)]
            out = out / chrom_norms
        }
        
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

    sample_types = rep(sample_classes, vapply(dt_list, nrow, numeric(1)))


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
                         aes(x = 0, y = 0, xend = PC1 * 2,
                             yend = PC2 * 2),
                         size = 1, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "black") + 
        annotate("text", x = pca_loadings$PC1 * 2.2,
                 y = pca_loadings$PC2 * 2.2, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab ("PC1") + ylab("PC2")


    return(p)

}


#' Plot substitution types
#' ### Two types of inputs
#'
#' 
#' 
#' @export


plot_sub_type_boxplots  <-  function(dt_list, sub, sample_classes, treatment) {

    plot_material = do.call( rbind,
                            lapply(names(dt_list),
                                   function(x) {
                                       df = dt_list[[x]][, sub, drop = FALSE];
                                       rownames(df) = paste(x, rownames(df),
                                                            sep = ":");
                                       return(df)
                                   })
                            )

    plot_material = data.frame(plot_material)
    colnames(plot_material) = "count"
    
    sample_types = rep(sample_classes, vapply(dt_list, nrow, numeric(1)))
    
    plot_material$smp = sample_types

    if (! missing(treatment)) {
        sample_treatment = rep(treatment, vapply(dt_list, nrow, numeric(1)))
        plot_material$treatment = sample_treatment
    }

    p = plot_count_boxplots(plot_material, treatment)
}


#' Plotting from a gg_df 
#' This part of code is wrapped in a separate function to
#' be used for both multinucleotide substitutions and SBSs
#'
#' 
#' 
#' @export


plot_count_boxplots <- function(plot_material, treatment) {
    if (! missing(treatment)){

        p = ggplot(plot_material, aes(x = smp, y = count, fill = treatment))
        
    } else {
        p = ggplot(plot_material, aes(x = smp, y = count)) 
    }

    p = p + geom_boxplot(outlier.shape = NA,
                         position = position_dodge(width = 0.8) )  +
        geom_jitter(position = position_dodge(width = 0.8))
    
    p = p + scale_color_brewer(palette = "Set1") + theme_bw() 
    
    return(p)
}


#' Plot all substitution plots
#' 
#'
#' 
#' 
#' @export


plot_sub_plots  <-  function(dt_list, sample_classes, treatment) {

    plot_list = lapply(
        MUT_TYPES, function(sub)
            plot_sub_type_boxplots(dt_list, sub, sample_classes, treatment ) +
            xlab("")
    )


    pout = ggarrange(plotlist=plot_list,nrow = 2, ncol = 3, labels = MUT_TYPES)
}



#' Plot counts
#' 
#'
#' 
#' 
#' @export


plot_counts  <-  function(mut_counts, sample_classes, treatment) {

    gg_df = data.frame(mut_counts)
    names(gg_df) = "counts"
    
    gg_df$classes = smp_classes

    if (! missing(treatment))  {
        gg_df$treatment = treatment

        p = ggplot(gg_df, aes(x= classes, y = counts, fill = treatment )) +
            geom_bar(stat = "identitfy", position = position_dodge(),
                     width = 0.8)
        
        
    } else {

        p = ggplot(gg_df, aes(x= classes, y = counts)) +
            geom_bar(stat = "identitfy", widht = 0.8)
    }
    p + xlab("") + ylab("Counts")
}




#' Returns number of multinucleotide substitutions for each sample
#' and each chromosome from a list of samples with elements as 
#' GRangesList with individual chromosome variants
#' 
#' 
#' @export


get_multinucl_profiles <- function(multisub_chr_list) {
    out = setNames(
        lapply(multisub_chr_list,
               function(sample_chrs)  {
                   lapply(sample_chrs, function(chr) {
                       vars = mcols(chr)[, c("REF", "ALT")]
                       
                  return(table(paste(vars$REF, vars$ALT, sep = ">") ) )
                   }
                  )
               }
               ), names(multisub_chr_list)
    )
    return(out)
}



#' Generates a count matrix with substitution types as columns
#' and samples-chrosomosomes as rows
#' 
#' 
#' 
#' @export



get_multisub_count_matrix <- function(multisub_chr_profile) {

    all_sub_types = c()

    for (sname in names(multisub_chr_profile))  {
        all_sub_types = c(all_sub_types,
                          unique(
                              do.call(c,
                                      lapply(multisub_chr_profile[[sname]], function(x)
                                          names(x)))) )
    }
    all_sub_types = unique(all_sub_types)
    
    all_sub_matrix = matrix(0, nrow = sum(sapply(multisub_chr_profile,length)),
                            ncol = length(all_sub_types))

    colnames(all_sub_matrix) = all_sub_types

    rownames(all_sub_matrix) = do.call(
        c, lapply( names(multisub_chr_profile),
                  function(x)
                      paste(x, names(multisub_chr_profile[[x]]), sep = ":")
                  ) )

    
    for (sname in names(multisub_chr_profile)) {
        for (chr in names(multisub_chr_profile[[sname]])) {
            s_c_subs = multisub_chr_profile[[sname]][[chr]] 
            all_sub_matrix[paste(sname, chr, sep = ":"), names(s_c_subs)] = s_c_subs
        }
    }
    
    return(all_sub_matrix)
}
