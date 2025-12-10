library(STdeconvolve)


run_stdecon <- function(counts_file, obs_file, results_dir_path, sample_name, k, ncores = 8, ngenes = 2000, alpha = 0.05) {

    # Read the counts and observation files
    cd <- read.csv(counts_file, sep = "\t", row.names=1)


    obs_df <- read.csv(obs_file, sep = "\t", row.names = 1)


    # # Fix the filename for subset_df
    subset_df <- read.csv(gsub("counts.txt", "subset_obs.txt", obs_file), sep = "\t", row.names = 1)
    # # Subset cd and obs_df based on the rownames of subset_df
    cd <- cd[rownames(cd) %in% rownames(subset_df), ]
    obs_df <- obs_df[rownames(obs_df) %in% rownames(subset_df), ]


    cd <- t(cd)
    cd <- round(cd)

    # Select highly variable genes and create a corpus
    pdf(file.path(results_dir_path, paste0("odg_stats_plot_lda_lda", sample_name, ".pdf")))
    counts <- cleanCounts(counts = cd, min.lib.size = 1, min.reads = 1, min.detected = 1, verbose = TRUE)
    corpus <- restrictCorpus(counts, nTopOD = ngenes, removeAbove = 1.0, removeBelow = 0.05, alpha = alpha, plot = TRUE, verbose = TRUE)
    writeLines(unlist(rownames(corpus)), file.path(results_dir_path, paste0("top_genes_stdecon_stdecon_", sample_name, ".csv")), sep = '\n')
    dev.off()
    
    mat <- t(as.matrix(corpus))
    zero_sum_rows <- apply(mat, 1, sum) == 0 
    mat <- mat[!zero_sum_rows, ]

    # # Train LDA model
    # print('Running fitLDA')
    ldas <- fitLDA(mat, Ks = k, perc.rare.thresh = 0.05, ncores = ncores, plot = FALSE, verbose = TRUE)

    # Save the results
    print(paste('Saving model with k =', k))
    optLDA <- optimalModel(models = ldas, opt = k)
    results <- getBetaTheta(optLDA, perc.filt = 0.025, betaScale = 1000)
    topics_per_spot <- results$theta
    genes_per_topic <- t(results$beta)

    # # Add prefix "lda_" to all column names
    colnames(topics_per_spot) <- paste0("stdecon_", colnames(topics_per_spot))
    colnames(genes_per_topic) <- paste0("stdecon_", colnames(genes_per_topic))
    
    # # Set row names of deconProp and deconGexp to be the index of obs_df
    # rownames(topics_per_spot) <- rownames(obs_df)  ##TODO Check 

    write.csv(topics_per_spot, file.path(results_dir_path, paste0("topics_per_spot_stdecon_", sample_name, ".csv")), row.names = TRUE)
    write.csv(genes_per_topic, file.path(results_dir_path, paste0("genes_per_topic_stdecon_", sample_name, ".csv")), row.names = TRUE)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
obs_file <- args[2]

results_dir_path <- args[3]
sample_name <- args[4]
k <- as.numeric(args[5])
ncores <- as.numeric(args[6])

print(args)
run_stdecon(counts_file, obs_file, results_dir_path, sample_name, k, ncores)
