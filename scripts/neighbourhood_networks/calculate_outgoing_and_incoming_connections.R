calculate_outgoing_and_incoming_connections <- function(graphobj){

    graph <- graphobj

    # Nodes and cluster info
    nodeDF <- data.frame(program = names(V(graph)), cluster = V(graph)$cluster)
    allclusters <- unique(nodeDF$cluster)
    # edge info 
    edgeDF <- data.frame(as_edgelist(graph))

    in_group_connections <- data.frame()
    out_group_connections <- data.frame()

    for(i in 1:nrow(nodeDF)){
        
        ref_program <- nodeDF$program[i]
        
        node_cluster <- nodeDF$cluster[which(nodeDF$program == ref_program)]

        cluster_members <- setdiff( nodeDF$program[which(nodeDF$cluster == node_cluster)], ref_program ) 

        in_group_tempDF <- data.frame(matrix(0, nrow = 1, ncol = length(allclusters), dimnames = list(ref_program, paste0("in_Cluster", allclusters))))
        out_group_tempDF <- data.frame(matrix(0, nrow = 1, ncol = length(allclusters), dimnames = list(ref_program, paste0("out_Cluster", allclusters)))) 
        
        edgeDF_in_outgoing <- edgeDF %>% filter(X1 == ref_program) %>% filter(X2 %in% cluster_members)
        edgeDF_in_incoming <- edgeDF %>% filter(X2 == ref_program) %>% filter(X1 %in% cluster_members)
        edgeDF_in <- rbind(edgeDF_in_outgoing,edgeDF_in_incoming) %>% unique()

        in_group_tempDF[,paste0("in_Cluster", node_cluster)] <- nrow(edgeDF_in)

        other_clusters <- setdiff(allclusters, node_cluster)
        for(j in 1:length(other_clusters)){
            
            other_cluster_members <- nodeDF$program[which(nodeDF$cluster == other_clusters[j])]
            
            edgeDF_out_outgoing <- edgeDF %>% filter(X1 == ref_program) %>% filter(X2 %in% other_cluster_members)
            edgeDF_out_incoming <- edgeDF %>% filter(X2 == ref_program) %>% filter(X1 %in% other_cluster_members)
            edgeDF_out <- rbind(edgeDF_out_outgoing, edgeDF_out_incoming) %>% unique()

            out_group_tempDF[,paste0("out_Cluster", other_clusters[j])] <- nrow(edgeDF_out)
        }

        in_group_connections <- rbind(in_group_connections, in_group_tempDF)
        out_group_connections <- rbind(out_group_connections, out_group_tempDF)
    }

    group_connections <- cbind(in_group_connections, out_group_connections)
    # return(group_connections)

    columnannot <- data.frame(V1 = colnames(group_connections))
    columnannot <- columnannot %>% separate(V1, c("connections", "cluster"), sep="_")
    rowannot <- nodeDF
     
    # cluster_cols <- pals::kelly()[3:(length(allclusters) + 2)]
    cluster_cols <- viridis(length(allclusters))
    names(cluster_cols) <- paste0("Cluster",allclusters)
    column_ha = HeatmapAnnotation(connections = columnannot$connections, cluster = columnannot$cluster, 
                                col = list(connections = c("in" = "black", "out" = "grey"), 
                                            cluster = cluster_cols))
                
    row_ha = rowAnnotation(cluster = paste0("Cluster", rowannot$cluster), 
                            col = list(cluster = cluster_cols))
    

    h1 <- Heatmap(as.matrix(group_connections), name = "n.edges", top_annotation = column_ha, left_annotation = row_ha, 
                column_split = columnannot$connections, row_split = paste0("Cluster", rowannot$cluster))
    h2 <- Heatmap(log2( as.matrix(group_connections) + 1), name = "log2(n.edges+1)", top_annotation = column_ha, left_annotation = row_ha, 
                column_split = columnannot$connections, row_split = paste0("Cluster", rowannot$cluster))
    
    data_list <- list("connections_matrix" = group_connections,
                    "h1" = h1, "h2" = h2)
    return(data_list)
}
