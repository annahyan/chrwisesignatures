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
#' @importFrom IRanges subsetByOverlaps
#' @importFrom philentropy H



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



mut_96_occurrences <- getFromNamespace("mut_96_occurrences", "MutationalPatterns")



#' Get mutational matrix in trinucleotide context for individual chromsomes
#' 
#'
#' 
#' 
#' @export




mut_matrix_chr  <-  function( chr_split_grangeslist, ref_genome, num_cores) {
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

    assembly_norm_factors = get_chrom_norm_factors(genome_assembly) 
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
            chrom_norms = assembly_norm_factors [rownames(smat)]
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


plot_coda_pca  <-  function(dt_list, sample_classes, point_size, arrow_length, method) {

    
    if (missing(point_size)) {
        point_size = 3
    }
        
    if (missing(arrow_length)) {
        arrow_length = 2
    }
    if (missing(method)) {
        method = "robust"
    }
    
    if (is.list(dt_list)) { 
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
    } else {
        plot_material = t(dt_list)
        sample_types = sample_classes
    }


    
    pca_list = pcaCoDa(plot_material, method = method)
    outliers = outCoDa(plot_material)
    print(outliers)


    ## Getting proporition of variance explained
    proportion_variance = (pca_list$eigenvalues ^ 2) / sum((pca_list$eigenvalues ^ 2))
    percent_variance = 100 * proportion_variance

    
    rob_pca_scores = data.frame(pca_list$scores)
    rob_pca_scores$col = sample_types

    rownames(rob_pca_scores) = rownames(plot_material)


    p = ggplot(rob_pca_scores, aes(x = Comp.1, y = Comp.2, color = col)) + 
        geom_point(size = point_size) +
        scale_color_brewer(palette = "Set1" ) # 
#        scale_color_manual(values = sample(CEMM_COLORS_ALL, 9 ) ) # +
    ## ggtitle(paste0("Rank:", nmf.rank))


    pca_loadings = as.data.frame(pca_list$loadings)
    rownames(pca_loadings) = colnames(plot_material)

    p = p + geom_segment(data = pca_loadings,
                         aes(x = 0, y = 0, xend = Comp.1 * arrow_length,
                             yend = Comp.2 * arrow_length),
                         size =0.5, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "gray50") + 
        annotate("text", x = pca_loadings$Comp.1 * arrow_length * 1.1,
                 y = pca_loadings$Comp.2 * arrow_length * 1.1, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab(paste0("PC1", " (",round(percent_variance[1], 2), "%)" ) )  +
        ylab(paste0("PC2", " (",round(percent_variance[2], 2), "%)") )
    
    return(p)
}

#' Plot CoDa PCA continuous coloring
#' 
#'
#' 
#' 
#' @export


plot_coda_pca_continuous  <-  function(dt_list, continuous, palette_name, point_size, arrow_length, method) {

    if (missing(point_size)) {
        point_size = 3
    }
    
    if (missing(arrow_length)) {
        arrow_length = 2
    }
    if (missing(method)) {
        method = "robust"
    }

    if (missing(palette_name)) {
        palette_name = "Spectral"
    }
    
    if (is.list(dt_list)) { 
        plot_material = do.call( rbind,
                                lapply(names(dt_list),
                                       function(x) {
                                           df = dt_list[[x]];
                                           rownames(df) = paste(x, rownames(df),
                                                                sep = ":");
                                           return(df)
                                       })
                                )

        ##        sample_types = rep(sample_classes, vapply(dt_list, nrow, numeric(1)))
    } else {
        plot_material = t(dt_list)
        ##        sample_types = sample_classes
    }


    
    pca_list = pcaCoDa(plot_material, method = method)
    outliers = outCoDa(plot_material)
    print(outliers)


    ## Getting proporition of variance explained
    proportion_variance = (pca_list$eigenvalues ^ 2) / sum((pca_list$eigenvalues ^ 2))
    percent_variance = 100 * proportion_variance

    
    rob_pca_scores = data.frame(pca_list$scores)
    rob_pca_scores$col = continuous

    rownames(rob_pca_scores) = rownames(plot_material)


    p = ggplot(rob_pca_scores, aes(x = Comp.1, y = Comp.2) ) + 
        geom_point(aes(colour = col), size = point_size) +
        scale_color_distiller(palette = palette_name ) # 
#        scale_color_manual(values = sample(CEMM_COLORS_ALL, 9 ) ) # +
    ## ggtitle(paste0("Rank:", nmf.rank))


    pca_loadings = as.data.frame(pca_list$loadings)
    rownames(pca_loadings) = colnames(plot_material)

    p = p + geom_segment(data = pca_loadings,
                         aes(x = 0, y = 0, xend = Comp.1 * arrow_length,
                             yend = Comp.2 * arrow_length),
                         size =0.5, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "gray50") + 
        annotate("text", x = pca_loadings$Comp.1 * arrow_length * 1.1,
                 y = pca_loadings$Comp.2 * arrow_length * 1.1, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab(paste0("PC1", " (",round(percent_variance[1], 2), "%)" ) )  +
        ylab(paste0("PC2", " (",round(percent_variance[2], 2), "%)") )
    
    return(p)
}


#' Plot classical PCA
#' 
#'
#' 
#' 
#' @export


plot_pca  <-  function(dt_list, sample_classes, point_size, arrow_length) {

        
    if (missing(point_size)) {
        point_size = 3
    }

        
    if (missing(arrow_length)) {
        arrow_length = 2
    }

    
    if (is.list(dt_list)) { 
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
    } else {
        plot_material = t(dt_list)
        sample_types = sample_classes
    }

    
    pca_list = prcomp(plot_material)

    percent_variance = 100 * summary(pca_list)$importance[2,]


    pca_scores = data.frame(pca_list$x)
    pca_scores$col = sample_types

    rownames(pca_scores) = rownames(plot_material)


    p = ggplot(pca_scores, aes(x = PC1, y = PC2, color = col)) + 
        geom_point(size = point_size) +
        scale_color_brewer(palette = "Set1" ) # 
#        scale_color_manual(values = sample(CEMM_COLORS_ALL, 9 ) ) # +
    ## ggtitle(paste0("Rank:", nmf.rank))


    pca_loadings = as.data.frame(pca_list$rotation)
    rownames(pca_loadings) = colnames(plot_material)

    p = p + geom_segment(data = pca_loadings,
                         aes(x = 0, y = 0, xend = PC1 * arrow_length,
                             yend = PC2 * arrow_length),
                         size = 1, 
                         arrow = arrow(length = unit(1/2, "picas") ),
                         color = "gray50") + 
        annotate("text", x = pca_loadings$PC1 * arrow_length * 1.1,
                 y = pca_loadings$PC2 * arrow_length * 1.1, 
                 label = rownames(pca_loadings)) +
        theme_bw() + theme(legend.title = element_blank()) +
        xlab(paste0("PC1", " (",round(percent_variance[1], 2), "%)" ) )  +
        ylab(paste0("PC2", " (",round(percent_variance[2], 2), "%)") )
    
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

    p = plot_count_boxplots(plot_material, "count", "smp", "treatment")
    return(p)
    
}


#' Plotting from a gg_df 
#' This part of code is wrapped in a separate function to
#' be used for both multinucleotide substitutions and SBSs
#'
#' 
#' 
#' @export


plot_count_boxplots <- function(plot_material, count_colname, sample, treatment) {
    if (! missing(treatment)){

        p = ggplot(plot_material, aes_string(x = sample,
                                             y = count_colname, fill = treatment))
        
    } else {
        p = ggplot(plot_material, aes(x = "smp", y = count_colname)) 
    }

    p = p + geom_boxplot(outlier.shape = NA,
                         position = position_dodge(width = 0.8) )  +
        geom_jitter(position = position_dodge(width = 0.8))
    
    p = p + # scale_color_brewer(palette = "Set1") +
        theme_bw() 
    
    return(p)
}


## #' Plot all substitution plots
## #' 
## #'
## #' 
## #' 
## #' @export


## plot_sub_plots  <-  function(dt_list, sample_classes, treatment) {

##     plot_list = lapply(
##         MUT_TYPES, function(sub)
##             plot_sub_type_boxplots(dt_list, sub, sample_classes, treatment ) +
##             xlab("")
##     )


##     pout = ggarrange(plotlist=plot_list,nrow = 2, ncol = 3, labels = MUT_TYPES)
## }



#' Plot counts
#' 
#'
#' 
#' 
#' @export


plot_counts  <-  function(mut_counts, sample_classes, treatment) {

    gg_df = data.frame(mut_counts)
    names(gg_df) = "counts"
    
    gg_df$classes = sample_classes

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



#' Split variants into genomic regions 
#' exons, introns, utr3, utr5, intergenic
#' 
#' 
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


#' MutationalPatterns mut_matrix added n_cores option
#' 
#' 
#' 
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


#' For genomic windows calculate Shannon entropy along with rainfall plots
#' 
#' 
#' 
#' 
#' @export

get_variant_entropy = function(variants, chr_length, window = 1e6, step = 500000) {

    step_pos = seq(1, chr_length - window + step, step)
    
    out = vapply(step_pos, function(pos) {
        relevant_starts = start(IRanges::subsetByOverlaps(ranges(variants),
                                                  IRanges(start = pos,
                                                         end = pos + window) ) )
        var_dists = log10(diff(relevant_starts))

        if (sum(var_dists) == 0) {
            return(0)
        } else {
            return(philentropy::H(var_dists / sum(var_dists)))
        }
    }, as.numeric(1))

    outdt = data.table(pos = step_pos, H = out)
    return(outdt)
}


#' Get regions with of FP cluster variants with given entropy_thresholds
#' and length of a region with high entropy.
#' 
#' 
#' 
#' 
#' @export


get_corrupt_regions = function(entropy_vals, threshold, window, step, expand) {

    if(missing(expand)) {
        expand = window
    }

    corrupt_regions = reduce(IRanges(entropy_vals[H > threshold, pos],
                                        width = step + 1))
    corrupt_regions = corrupt_regions[width(corrupt_regions) > window]


    out = resize(corrupt_regions, width = width(corrupt_regions) + expand,
                 fix = "start")

    return(out)
    
}



#' CoDa correlation of all pairs with p.values 
#' 
#' Slightly modified from the function from robCompositions
#' 
#' @export


corCoDa = function(x, p.val = 0.05, ...) {
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    if (any(x[!is.na(x)] <= 0)) 
        stop("all elements of x must be greater than 0")
    if (ncol(x) <= 2) 
        stop("calculation of average symmetric coordinates not possible")
    balZav <- function(x) {
        D <- ncol(x)
        Z.av <- matrix(NA, ncol = 2, nrow = nrow(x))
        p1 <- sqrt(D - 1 + sqrt(D * (D - 2)))/sqrt(2 * D)
        if (D == 3) {
            p2 <- x[, 3]
        }
        else {
            p2 <- apply(x[, 3:D], 1, prod)
        }
        p3 <- (sqrt(D - 2) + sqrt(D))/(sqrt(D - 2) * (D - 1 + 
            sqrt(D * (D - 2))))
        p4 <- 1/(D - 1 + sqrt(D * (D - 2)))
        Z.av[, 1] <- p1 * (log(x[, 1]/(x[, 2]^p4 * p2^p3)))
        Z.av[, 2] <- p1 * (log(x[, 2]/(x[, 1]^p4 * p2^p3)))
        return(Z.av)
    }
    ind <- c(1:ncol(x))
    corZav <- matrix(NA, ncol(x), ncol(x))
    corPvals <- matrix(NA, ncol(x), ncol(x))
    for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {
            balZavout = balZav(x[, c(i, j, ind[-c(i, j)])])
            corPvals[i, j]  <-  cor.test(balZavout[,1], balZavout[, 2], ...)$p.value

            corZav[i, j] <- cor(balZavout, ...)[1, 2]
        }
    }
    corZav[corPvals > p.val ] = 0
    corZav[lower.tri(corZav)] <- t(corZav)[lower.tri(corZav)]
    diag(corZav) <- 1
    return(corZav)
}




