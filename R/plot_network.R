plot_network <- function(adjacency_matrix, min_threshold = 0.2, layout) {

    if (is.missing(layout) ) {
        layout = "kk"
    }


    graph.input = adjacency_matrix
    graph.input[ abs(graph.input) < min_threshold ] = 0
    graph.input[graph.input == 1] = 0
    graph.input[abs(graph.input) >= min_threshold ] = 1

    ## adjacency_matrix[ abs(adjacency_matrix) < 0.2 ]

    sig.adjacency = graph.adjacency(graph.input)
    
    graph <- graph_from_data_frame(get.edgelist(sig.adjacency) )


    V(graph)$label = V(graph)$name
    E(graph)$intensity = sapply( 1:length(E(graph)), 
                                function(x) {
                                    gends = ends(graph, x)
                                    adjacency_matrix[gends[1], gends[2] ] 
                                }
                                )
    E(graph)$edgewidth = abs(E(graph)$intensity) 

    P = ggraph(graph, layout = layout) + 
        geom_edge_link(aes(color = intensity, width = edgewidth)) + 
        scale_edge_color_gradient2(low = rgb(0, 140, 160, maxColorValue = 255), 
                                   high = rgb(210, 50, 60, maxColorValue = 255)) + 
        scale_edge_width(range = c(0.8, 1.3), guide = FALSE) +
        geom_node_point() + 
        expand_limits(x = c(0, 4))+
        geom_node_label(aes(label = label)) +
        theme_void() +
        theme(legend.position = "right", legend.title =element_blank() )
    return(P)
}
