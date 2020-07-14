#' Lineplot variant counts across chromosomes and samples
#'
#' @details
#' Plot chromosome-wise counts with samples as facets, colored by conditions,
#' clones are grouped.
#'
#' @param chr_wise_counts count matrix, result from count_chr_variants
#' @param KOs sample KO's by which plot facets are defined
#' @param treatment treatment conditions by which coloring is defined
#' @param clones replicates 
#'
#' @import ggplot2
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
