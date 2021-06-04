#' Interactive signature network from adjacency-like matrix
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
#' @param layout Layout passed to visNetwork. Default: layout_nicely.
#'
#' 
#' @import igraph
#' @import ggraph
#' @import visNetwork
#' 
#' @export


plot_interactive_network <- function(adjacency_matrix,
                                     min_threshold = 0.2, binary_matrix = FALSE, layout) {


    if (missing(layout) ) {
        layout_custom = "layout_nicely"
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

    data = toVisNetworkData(graph)

    nodes = data$nodes

    nodes$font.size = 25
    
    edges = data$edges

    
    edges$intensity = sapply( 1:length(E(graph)),
                                function(x) {
                                    gends = ends(graph, x)
                                    adjacency_matrix[gends[1], gends[2] ]
                                } )

    edges$width = abs(edges$intensity)

    
    if (binary_matrix | nlevels(factor(E(graph)$intensity)) == 1 ) {
        edges$color = "black"
    } else {

        col.edges = c( rgb(0, 140, 160, maxColorValue = 255),
                      rgb(210, 50, 60, maxColorValue = 255) )
        
        edges$color =  col.edges[ (edges$intensity > 0) + 1 ] 
            }

    visNetwork(nodes, edges) %>% visIgraphLayout(layout = layout_custom) %>%
        visLegend()
}

