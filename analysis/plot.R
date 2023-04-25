setwd("/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/plot")
source("plot_utils.R")

# File path where GCN output stored
path_input <- "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/"
path_tf_and_reqdgenes <- "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/"
TEST_ID <- "TEST4"

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
    list_predicted[[TYPE]] <- read_graph_prediction(TYPE, 'TEST3', cutoff = 0.8)
}
names(list_predicted)[3]= 'GAT_pseudo'
list_predicted[['GAT']] <- read_graph_prediction(TYPE, 'TEST4', cutoff = 0.8)
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
plot_f1(results_dorothea, TEST_ID, "Dorothea")
plot_f1(results_gtrd, TEST_ID, "GTRD")

# plot number of edges
plot_n_edges(networks)
