library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(ComplexHeatmap)
library(parallel)
library(future)
library(future.apply)
library(doParallel) 
library(osfr)
library(circlize)
library(RColorBrewer)
library(openxlsx)
library(scico)
library(stats)
library(rstatix)
library(cowplot)


cc_interaction_networks <- function(usage_matrix, metadata, metadata_column, cond1, cond2, save_folder, save_subfolder, n_bins, GLOBAL_THRESHOLD, tumorTME = FALSE, fileprefix = NULL, edge_threshold, pattern_regex, replace_with, plot_pval_heatmapconnections_heatmap){

    spotmetadata = metadata
    usage_matrix = usage_matrix[rownames(spotmetadata),]

    SAVE_FOLDER = save_folder

    CC_INTERACTION_ENRICHMENT = TRUE
    prefix = paste0(fileprefix, edge_threshold, "_")
    
    ## CC_INTERACTION_ENRICHMENT

    # rm(list = setdiff(setdiff(ls(), lsf.str()), program_var_names))
    save_path <- paste(SAVE_FOLDER, "/cell_cell_interactions/", sep = "")
    dir.create(save_path)
    dir.create(paste0(save_path, "samples_split/"))

    cutoff <- as.data.frame(usage_matrix)
    # colnames(cutoff) <- program_annotations
    threshold <- GLOBAL_THRESHOLD

    if(file.exists(paste(save_path,  "uthresh", threshold, "_", fileprefix,  "_Cell_Cell_Interaction_Enrichment.tsv", sep = ""))){
        stats_df <- read.csv(paste(save_path,  "uthresh", threshold, "_", fileprefix, "_Cell_Cell_Interaction_Enrichment.tsv", sep = ""), sep = "\t")
    } else{
        # Set up parallel backend
        n_cores <- parallelly::availableCores() - 1
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)

        cutoff_rows <- rownames(cutoff)
        programs <- colnames(cutoff)

        stats_df <- foreach(p1 = programs, .combine = rbind, .packages = c("stats", "rstatix")) %:%
        foreach(p2 = programs, .combine = rbind) %dopar% {

            above_threshold <- which(cutoff[, p1] > threshold)
            p1_bins_above_threshold <- cutoff_rows[above_threshold]
            
            p1_cond1_bins <- p1_bins_above_threshold

            cond1_vals <- cutoff[p1_cond1_bins, p2, drop = FALSE]

            contingency_table <- matrix(c(sum(cond1_vals > threshold, na.rm = TRUE),
                                        sum(cond1_vals <= threshold, na.rm = TRUE)),
                                        nrow = 1,
                                        byrow = FALSE,
                                        dimnames = list(c(cond1), c("p2_pos", "p2_neg")))

            #ft_df <- pairwise_fisher_test(contingency_table, p.adjust.method = "fdr")
            
            ft_df <- data.frame(group1 = cond1, 
                                n = sum(contingency_table))
            ft_df$program_one <- p1
            ft_df$program_two <- p2
            ft_df[[paste0(cond1, "_P2pos")]] <- contingency_table[cond1, "p2_pos"]
            ft_df[[paste0(cond1, "_P2neg")]] <- contingency_table[cond1, "p2_neg"]

            ft_df
        }
        stopCluster(cl)

        write.table(stats_df, paste(save_path,  "uthresh", threshold, "_", fileprefix, "_Cell_Cell_Interaction_Enrichment.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
    }

    # 4.1 Network plot of the data ----
    ## 4.1.1 Modify the edge values ----
    # Add pseudocount and compute sample_idal probabilities
    add_pseudocount <- function(pos, neg) (pos + 1) / (pos + 1 + neg + 1)
    stats_df[[paste0(cond1, "_val")]]  <- add_pseudocount(stats_df[[paste0(cond1, "_P2pos")]], stats_df[[paste0(cond1, "_P2neg")]])

    stats_df <- stats_df %>% mutate(gt01 = as.integer(abs(stats_df[[paste0(cond1, "_val")]]) > 0.099),
                                    gt02 = as.integer(abs(stats_df[[paste0(cond1, "_val")]]) > 0.199),
                                    gt03 = as.integer(abs(stats_df[[paste0(cond1, "_val")]]) > 0.299),
                                    gt04 = as.integer(abs(stats_df[[paste0(cond1, "_val")]]) > 0.399),
                                    gt05 = as.integer(abs(stats_df[[paste0(cond1, "_val")]]) > 0.499))

    cutoff$sample <- cond1 

    node_summary <- cutoff %>%
        mutate(sample_id = cond1) %>%
        group_by(sample_id) %>%
        mutate(total_spots = n()) %>%
        group_by(sample, sample_id, total_spots) %>%
        summarise(across(starts_with("ot"), ~ sum(. > threshold), .names = "{col}"), .groups = "drop") %>%
        pivot_longer(cols = starts_with("ot"), names_to = "program", values_to = "cps") %>%
        group_by(program) %>%
        mutate(num_samples_cps_gt_100 = sum(cps >= 100)) %>%
        group_by(sample_id, program) %>%
        mutate(cpp = sum(cps)) %>%
        mutate(proportion = cpp / total_spots) %>%
        group_by(program) %>%
        ungroup() %>%
        select(sample_id, total_spots, program, num_samples_cps_gt_100, cpp, proportion) %>%
        distinct()

    for (program in unique(node_summary$program)){
        cond1_proportion <- node_summary$proportion[node_summary$program == program & node_summary$sample_id == cond1]
    }



    ## Generate the network ----
    node_attrs <- unique(node_summary[, c("program", "num_samples_cps_gt_100")])

    # First sample specific network
    for (sample in paste0(cond1,"_val")){

        edges <- data.frame(from = stats_df$program_one,
                          to = stats_df$program_two,
                          weight = abs(stats_df[[sample]]),
                          weight_col = stats_df[[sample]]) #,
                          #sig_edge = stats_df$p)
    
        # 1. Create the igraph object (directed) ----
        graph <- graph_from_data_frame(d = edges, vertices = node_attrs, directed = TRUE)
        # stats_df$edge_threshold <- F; stats_df$edge_threshold[stats_df[[edge_threshold]] == 1 & stats_df$n >= n_bins & stats_df$p <= 0.05] <- T
        stats_df$edge_threshold <- F; stats_df$edge_threshold[stats_df[[edge_threshold]] == 1 & stats_df$n >= n_bins ] <- T
        edge_indicies <- which(stats_df$edge_threshold == F)
        same_nodes <- which(stats_df$program_one == stats_df$program_two)
        graph <- delete_edges(graph, unique(c(edge_indicies,same_nodes)))

        # 2. Perform community detection using Infomap ----
        clusters <- cluster_infomap(graph)
        V(graph)$cluster <- clusters$membership
        V(graph)$color <- as.factor(clusters$membership)
        # V(graph)$diff_sign <- ifelse(V(graph)$fold_change >= 0, "Positive", "Negative")

        # 3. Make a graph and layout without the isolated nodes ----
        isolated_nodes <- which(degree(graph) == 0)
        graph_filtered <- delete_vertices(graph, isolated_nodes)

    

        # 6. Make the layouts ----
        layout_df <- create_layout(graph, layout = "graphopt")
        layout_filtered <- create_layout(graph_filtered, layout = "graphopt")
      
        ## Plotting ----
        base_plot <- function(graph_obj, layout_obj, node_color){
            hawaii_colors <- scico(n = 256, palette = "hawaii", direction = -1)

            ggraph(layout_obj) + theme_void() +
            geom_edge_fan(aes(edge_alpha = weight, edge_width = weight, edge_colour = weight_col),
                            arrow = arrow(length = unit(2, 'mm')), start_cap = circle(3, "mm"), end_cap = circle(3, "mm"),
                            show.legend = TRUE) +
            geom_node_point(aes(size = num_samples_cps_gt_100, color = {{ node_color }})) +
            geom_node_text(aes(label = name), repel = TRUE, size = 5, bg.colour = "white", bg.r = 0.1) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
                    plot.subtitle = element_text(hjust = 0.5, size = 14),
                    legend.title = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    plot.margin = margin(1, 1, 1, 1, "cm")) +
            scale_edge_colour_gradient(low = "black", high = "magenta",
                                        name = paste0(sample, " [Edge Color]")) +
            scale_edge_width_continuous(range = c(0.3, 3),
                                        name = paste0("| ", sample, " [Edge Color] |")) + # "-log10(P-Value) Pairwise\nFisher Exact Test [Edge Width]") +
            scale_edge_alpha_continuous(range = c(0.2, 1),
                                        name = paste0("| ", sample, " [Edge Color] |")) +
            scale_size_continuous(name = paste0("# of Samples With >= ", n_bins," Bins\nAbove 0.05 Usage [Node Size]"))
        }

        ### 1. Plot the whole network ----
        network_plot <- base_plot(graph, layout_df, V(graph)$fold_change) +
            scale_color_scico(palette = "berlin", direction = -1, limits = c(-5, 5),
                            name = "Symmetric Fold Change\n[Node Color Capped (-5, 5)]") +
            labs(title = paste0("Program-Program Interaction Network: ", sample),
                subtitle = paste("Interactions with |", sample, "| >", as.numeric(gsub("gt","",edge_threshold))/10, "& Bins >= ", n_bins," & P-Value <= 0.05
                                \n", cond1," Prevalent Programs Are Positive"))
        network_plot_wo_isolated_nodes <- base_plot(graph_filtered, layout_filtered, V(graph_filtered)$fold_change) +
            scale_color_scico(palette = "berlin", direction = -1, limits = c(-5, 5),
                            name = "Symmetric Fold Change\n[Node Color Capped (-5, 5)]") +
            labs(title = paste0("Program-Program Interaction Network: ", sample),
                subtitle = paste("Interactions with |", sample, "| >", as.numeric(gsub("gt","",edge_threshold))/10, "& Bins >= ", n_bins," & P-Value <= 0.05
                                \n", cond1," Prevalent Programs Are Positive
                                \nIsolated Nodes Removed"))
        if(max(table(V(graph)$cluster)) != 1){                       
            network_plot_wo_isolated_nodes_cluster_split_col <- network_plot_wo_isolated_nodes + facet_graph(~ cluster)
        }

        ### 2. Plot the whole network (clustered) ----
        network_plot_cluster <- base_plot(graph, layout_df, as.factor(V(graph)$cluster)) +
            scale_color_viridis_d(name = "Cluster") +
            labs(title = "Network Colored by Cluster",
                subtitle = paste("Interactions with |", sample, "| >", as.numeric(gsub("gt","",edge_threshold))/10, "& Bins >= ", n_bins," & P-Value <= 0.05
                                \n", cond1," Prevalent Programs Are Positive"))
        network_plot_wo_isolated_nodes_cluster <- base_plot(graph_filtered, layout_filtered, as.factor(V(graph_filtered)$cluster)) +
            scale_color_viridis_d(name = "Cluster") +
            labs(title = "Network Colored by Cluster",
                subtitle = paste("Interactions with |", sample, "| >", as.numeric(gsub("gt","",edge_threshold))/10, "& Bins >= ", n_bins," & P-Value <= 0.05
                                \n", cond1," Prevalent Programs Are Positive
                                \nIsolated Nodes Removed"))
        if(max(table(V(graph)$cluster)) != 1){
            network_plot_wo_isolated_nodes_cluster_split <- network_plot_wo_isolated_nodes_cluster + facet_graph(~ cluster)
        }


        ### Print the plot ----
        pdf(file.path(save_path, paste0("samples_split/", "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "nbins", n_bins,"_", sample, "_Program_Program_Interaction_Network_Edge_Cutoff.pdf")), width = 20, height = 15) # was w = 12, h = 10
        print(network_plot)
        print(network_plot_wo_isolated_nodes)
        print(network_plot_cluster)
        print(network_plot_wo_isolated_nodes_cluster)
        if(max(table(V(graph)$cluster)) != 1){
            print(network_plot_wo_isolated_nodes_cluster_split)
            print(network_plot_wo_isolated_nodes_cluster_split_col)
        }
        dev.off()

        saveRDS(graph, file.path(save_path, paste0("samples_split/", "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "nbins", n_bins,"_", sample, "_Program_Program_Interaction_Network_Edge_Cutoff__graphobject.rds")))
        saveRDS(layout_df, file.path(save_path, paste0("samples_split/", "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "nbins", n_bins,"_", sample, "_Program_Program_Interaction_Network_Edge_Cutoff__layoutobject.rds")))
        saveRDS(layout_filtered, file.path(save_path, paste0("samples_split/", "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "nbins", n_bins,"_", sample, "_Program_Program_Interaction_Network_Edge_Cutoff__filteredlayoutobject.rds")))
    }
    
    if(connections_heatmap){
        connections_data <- calculate_outgoing_and_incoming_connections(graphobj = graph)
        connec_matrix <- connections_data[["connections_matrix"]]
        connec_matrix <- data.frame(connec_matrix) %>% add_column(program = rownames(connec_matrix), .before = colnames(connec_matrix)[1]) 
        write.table(connec_matrix, file.path(save_path, paste0("samples_split/", "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "nbins", n_bins,"_", sample, "_Program_Program_Interaction_Network_Edge_Cutoff__connection_matrix.txt")), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

        pdf(paste0("samples_split/", "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "nbins", n_bins,"_", sample, "_Program_Program_Interaction_Network_Edge_Cutoff__connection_matrix.pdf"), height = 12, width = 8)
        draw(connections_data[["h1"]])
        draw(connections_data[["h2"]])
        dev.off()
    }


    if(plot_pval_heatmap){
        # Heatmap of the p-values ----
        # Columns are program one, and rows and program two
        stats_df_spread <- stats_df %>%
            select(program_one, program_two, p) %>%
            pivot_wider(names_from = program_one, values_from = p, values_fill = NA) %>%
            arrange(program_two)
        stats_df_spread <- as.data.frame(stats_df_spread)
        rownames(stats_df_spread) <- stats_df_spread$program_two
        stats_df_spread <- stats_df_spread[,-1]
        stats_df_spread <- as.matrix(stats_df_spread)
        stats_df_spread_cap <- -log10(stats_df_spread)
        stats_df_spread_cap[is.infinite(stats_df_spread_cap)] <- 0
        stats_df_spread_cap[stats_df_spread_cap >= -log10(0.0001)] <- 4

        htp_name <- paste("Fisher Exact P-Value ",cond1," vs ",cond2," [Post-Hoc Pairwise Fisher Exact P-Value]", sep = "")
        cols <- colorRamp2(breaks = range(stats_df_spread, na.rm = TRUE), hcl_palette = "Blues")
        cols_cap <- colorRamp2(breaks = range(stats_df_spread_cap, na.rm = TRUE), hcl_palette = "Reds", reverse = T)

        htp1 <- Heatmap(stats_df_spread, cluster_columns = FALSE, cluster_rows = FALSE,
                    show_row_names = TRUE, show_column_names = TRUE,
                    na_col = "white", col = cols, column_title = "Program One", row_title = "Program Two")
        htp2 <- Heatmap(stats_df_spread, cluster_columns = TRUE, cluster_rows = FALSE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols, column_title = "Program One", row_title = "Program Two")
        htp3 <- Heatmap(stats_df_spread, cluster_columns = FALSE, cluster_rows = TRUE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols, column_title = "Program One", row_title = "Program Two")
        htp4 <- Heatmap(stats_df_spread, cluster_columns = TRUE, cluster_rows = TRUE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols, column_title = "Program One", row_title = "Program Two")

        htp5 <- Heatmap(stats_df_spread_cap, cluster_columns = FALSE, cluster_rows = FALSE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols_cap, column_title = "Program One", row_title = "Program Two")
        htp6 <- Heatmap(stats_df_spread_cap, cluster_columns = TRUE, cluster_rows = FALSE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols_cap, column_title = "Program One", row_title = "Program Two")
        htp7 <- Heatmap(stats_df_spread_cap, cluster_columns = FALSE, cluster_rows = TRUE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols_cap, column_title = "Program One", row_title = "Program Two")
        htp8 <- Heatmap(stats_df_spread_cap, cluster_columns = TRUE, cluster_rows = TRUE,
                        show_row_names = TRUE, show_column_names = TRUE,
                        na_col = "white", col = cols_cap, column_title = "Program One", row_title = "Program Two")
        
        p11 <- grid.grabExpr(draw(htp1, padding = unit(c(3, 3, 3, 3), "cm")))
        p12 <- grid.grabExpr(draw(htp2, padding = unit(c(3, 3, 3, 3), "cm")))
        p13 <- grid.grabExpr(draw(htp3, padding = unit(c(3, 3, 3, 3), "cm")))
        p14 <- grid.grabExpr(draw(htp4, padding = unit(c(3, 3, 3, 3), "cm")))
        p21 <- grid.grabExpr(draw(htp5, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "-log10(p-value) [Capped at -log10(0.0001)]"))
        p22 <- grid.grabExpr(draw(htp6, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "-log10(p-value) [Capped at -log10(0.0001)]"))
        p23 <- grid.grabExpr(draw(htp7, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "-log10(p-value) [Capped at -log10(0.0001)]"))
        p24 <- grid.grabExpr(draw(htp8, padding = unit(c(3, 3, 3, 3), "cm"), column_title = "-log10(p-value) [Capped at -log10(0.0001)]"))

        pdf(paste(save_path, "usagethresh", GLOBAL_THRESHOLD, "_", prefix, "_Cell_Cell_Interaction_Enrichment_Nominal_P_Value.pdf", sep = ""), width = 80, height = 20)
        print(plot_grid(p11, p12, p13, p14, ncol = 4))
        print(plot_grid(p21, p22, p23, p24, ncol = 4))
        dev.off()
    }
}