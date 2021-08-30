#' Signature network from adjacency-like matrix
#'
#' @details 
#' Network representation of signature interaction matrix. The edges are colored
#' based on interaction sign and value.
#'
#' @param adjacency_matrix A sparse matrix with values values indicating the
#' strength and direction of interaction. Correlations, adjacencies, etc.
#'
#' @param min_threshold All values with less absolute value than min_threshold
#' will be set to 0. Default: 0.2
#'
#' @param binary_matrix If binary, the links are constructed from values 1.
#' Otherwise the matrix is considered a correlation matrix and all 0 values are
#' set to 0 before network construction. Default: FALSE.
#' 
#' @param layout Layout passed to ggraph. Default: kk.
#'
#' @param expand_limits Argument to control the expansion across x axis.
#' Default: c(0, 1)
#' 
#' @import igraph
#' @import ggraph
#' 
#' @export


plot_network <- function(adjacency_matrix,
                         min_threshold = 0.2,
                         binary_matrix = FALSE,
                         layout,
                         expand_limits) {

    if (missing(layout) ) {
        layout = "kk"
    }
    
    if (missing(expand_limits)) {
        expand_limits = c(0, 1)
    }
    
    if (all( abs(adjacency_matrix[upper.tri(adjacency_matrix) ] )  <=  min_threshold, na.rm = TRUE ) ) {
        return(ggplot() + theme_void ()) # return empty ggplot
    }

    graph.input = adjacency_matrix
    graph.input[ abs(graph.input) < min_threshold ] = 0

    ## if( ! binary_matrix ) {
    ##     graph.input[graph.input == 1] = 0
    ## }

    diag(graph.input) = 0
    
    graph.input[abs(graph.input) >= min_threshold ] = 1

    ## adjacency_matrix[ abs(adjacency_matrix) < 0.2 ]

    sig.adjacency = graph.adjacency(graph.input)

    graph <- graph_from_data_frame(get.edgelist(sig.adjacency) )

    V(graph)$label = V(graph)$name
    E(graph)$intensity = sapply( 1:length(E(graph)),
                                function(x) {
                                    gends = ends(graph, x)
                                    adjacency_matrix[gends[1], gends[2] ] 
                                } )
    
    E(graph)$edgewidth = abs(E(graph)$intensity) 

    P = ggraph(graph, layout = layout)
    if ( binary_matrix | nlevels(factor(E(graph)$intensity)) == 1  ) {
        if ( unique(E(graph)$intensity) == 1 ) {    
            col = rgb(210, 50, 60, maxColorValue = 255) # 1
        } else {
            
            col = rgb(0, 140, 160, maxColorValue = 255) # -1
        }
        P = P +  geom_edge_link(color = col)
    } else {
        P = P + geom_edge_link(aes(color = intensity, width = edgewidth)) + 
            scale_edge_color_gradient2(low = rgb(0, 140, 160, maxColorValue = 255), 
                                       high = rgb(210, 50, 60, maxColorValue = 255)) + 
            scale_edge_width(range = c(0.8, 1.3), guide = "none") +
            expand_limits( x = expand_limits)
        
            }
    P = P + geom_node_point() + 
        geom_node_label(aes(label = label)) +
        theme_void() +
        theme(legend.position = "right", legend.title =element_blank() )
    return(P)
    
}
