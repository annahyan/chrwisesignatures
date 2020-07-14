#' Classical PCA without CoDa transformation for comparison
#'
#' @details 
#' This function is provided solely as a comparison tool agains
#' \code{\link{plot_coda_pca}} \code{\link{plot_coda_pca_continuous}}.
#' 
#' @param dt_list Either list of mutational profiles or a data.frame with
#' compositional parts(signature intensities, mutation counts) as columns and
#' samples as rows.
#'
#' @param sample_classes A vector of sample classes to be color coded. The
#' length of this parameter has to be equal to the length of dt_list if list is
#' provided, or number of rows if data.frame provided.
#'
#' @param point_size Point size
#'
#' @param arrow_length the length of biplot arrows
#'
#' @importFrom robCompositions pcaCoDa 
#'
#' @seealso \code{\link{plot_coda_pca}}, \code{\link{plot_coda_pca_continuous}}.
#' 
#' @export


plot_pca  <-  function(dt_list, sample_classes, point_size, arrow_length) {

        
    if (missing(point_size)) {
        point_size = 3
    }

        
    if (missing(arrow_length)) {
        arrow_length = 2
    }

    
    if (inherits(dt_list, "list")) { 
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
