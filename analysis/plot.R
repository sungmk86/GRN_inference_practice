setwd("/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/analysis")
source("plot_utils.R")

# File path where GCN output stored
path_input <- "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/"
path_tf_and_reqdgenes <- "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/"

# get all genes
input_expr <- file.path(path_input, "Bayesian_DE/iterative_test/iterative_test/TF_experiment_expression_matrix.gz")
df_expr_full <- read.delim(input_expr, row.names = 1)
genes_all <- df_expr_full$id

# get all TFs
external_tfs <- grep("tf", genes_all, value = T)
tfs <- get_tfs(path_tf_and_reqdgenes)
tfs_all <- unique(c(external_tfs, tfs))


# Load graph predictions
list_predicted <- list()
for (TYPE in c("GCN", "graphSAGE", "GAT")) {
    if (TYPE == "GAT") {
        TEST_ID <- "TEST3"
        lst_rng_seed <- c("", "_rng111_neg2.0", "_rng123_neg2.0", "_rng1234_neg2.0", "_rng111_neg3.0", "_rng123_neg3.0", "_rng1234_neg3.0", "_rng111_neg5.0", "_rng123_neg5.0", "_rng1234_neg5.0", "_rng111_neg3.0_super", "_rng123_neg3.0_super", "_rng1234_neg3.0_super", "_rng111_neg3.0_linear", "_rng123_neg3.0_linear", "_rng1234_neg3.0_linear", "_rng111_neg3.0_suplinear", "_rng123_neg3.0_suplinear", "_rng1234_neg3.0_suplinear")

        for (rng in lst_rng_seed) {
            ID <- paste0(TEST_ID, rng)
            list_predicted[[paste0(TEST_ID,'_', TYPE, rng)]] <- read_graph_prediction(TYPE, ID)
        }
        TEST_ID <- "TEST4"
        lst_rng_seed <- c("_rng111_neg3.0_linear", "_rng123_neg3.0_linear", "_rng1234_neg3.0_linear")
        for (rng in lst_rng_seed) {
            ID <- paste0(TEST_ID, rng)
            list_predicted[[paste0(TEST_ID,'_', TYPE, rng)]] <- read_graph_prediction(TYPE, ID)
        }
        TEST_ID <- "TEST5"
        for (rng in lst_rng_seed) {
            ID <- paste0(TEST_ID, rng)
            list_predicted[[paste0(TEST_ID,'_', TYPE, rng)]] <- read_graph_prediction(TYPE, ID, .9)
        }
    } else {
        list_predicted[[TYPE]] <- read_graph_prediction(TYPE, "TEST3")
    }
}
networks <- load_all_networks(list_predicted, tfs_all)

df_all_possible_edges <- expand.grid(tfs_all, genes_all)
all_g <- igraph::graph_from_edgelist(as.matrix(df_all_possible_edges[, 1:2]))

# Reference
reference <- attr(E(simplify(all_g)), "vnames")
refers <- function(v) factor(ifelse(reference %in% attr(E(v), "vnames"), "Edge", "Not"), levels = c("Edge", "Not"))
ref_dorothea <- refers(networks$Dorothea$ig)
ref_GTRD <- refers(networks$GTRD$ig)

results_dorothea <- lapply(networks, function(x) {
    caret::confusionMatrix(refers(x$ig), ref_dorothea, mode = "everything")
})
results_gtrd <- lapply(networks, function(x) {
    caret::confusionMatrix(refers(x$ig), ref_GTRD, mode = "everything")
})

# Represent metrics
TEST_ID = 'TEST4'
plot_f1(results_dorothea, paste0(TEST_ID, "_rng_boxplot"), "Dorothea")
plot_f1(results_gtrd, paste0(TEST_ID, "_rng_boxplot"), "GTRD")

plot_score_distribution(dir_network = paste0("../GAT/data/"), TEST_ID)

plot_scores_comp(dir_network = paste0("../GAT/data/"), TEST_ID, lst_rng_seed, type = "eigen_centrality", tfs_all)
plot_edge_scores_comp(dir_network = paste0("../GAT/data/"), TEST_ID, lst_rng_seed)


# plot number of edges
# plot_n_edges(networks)
