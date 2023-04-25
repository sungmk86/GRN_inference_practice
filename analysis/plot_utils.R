library(igraph)
library(ggplot2)
library(data.table)
##########################

edge_type <- function(nt, tfs, index = c("source", "target")) {
  r <- paste(ifelse(nt$source %in% tfs, "TF", "target"),
    ifelse(nt$target %in% tfs, "TF", "target"),
    sep = "-"
  )
  return(r)
}

# Wrapper to read network, assign edge type, filter, and create and igraph object
read_network <- function(path_network, tfs_all, directed = FALSE, filter = c("TF-TF", "TF-target")) {
  if (class(path_network) == "data.frame") {
    net <- path_network
  } else {
    net <- fread(path_network, sep = "\t", header = F)
  }
  colnames(net) <- c("source", "target", "score", "is_directional")[1:ncol(net)]
  net <- net[!net$source %in% c("TF", "Tf", "tf"), ]
  net$type <- edge_type(net, tfs_all)
  # Filter desired edges
  net <- net[net$type %in% filter, ]
  # igraph object
  ig <- igraph::graph_from_edgelist(as.matrix(net[, c("source", "target")]), directed = directed)
  ig <- igraph::simplify(ig)
  return(list(edges = net, ig = ig))
}

load_all_networks <- function(list_predicted, tfs_all) {
  path_networks <- list(
    GCN = list_predicted[["GCN"]],
    graphSAGE = list_predicted[["graphSAGE"]],
    GAT_pseudo = list_predicted[["GAT_pseudo"]],
    GAT = list_predicted[["GAT"]],
    eBIC = "/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/BIC/edge_strength_filt.txt",
    STRING = "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/BIC/data/networks_anonymize.txt",
    GTRD = "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/gtrd_anonymized.txt",
    Dorothea = "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/dorothea_anonymized.txt"
  )

  networks <- sapply(path_networks, function(x) read_network(x, tfs_all, directed = T), simplify = F)
  return(networks)
}


plot_f1 <- function(results_f1, TEST_ID, type) {
  fig_name <- paste0("data/", TEST_ID, "_", type, ".pdf")
  accuracies <- sapply(results_f1, function(x) x$byClass)
  y <- melt(data.frame(metric = rownames(accuracies), accuracies))

  p <- ggplot(
    y[y$metric %in% c("F1", "Recall", "Precision") & y$variable != type, ],
    aes(x = variable, y = value * 100, fill = variable)
  ) +
    ylab("(%)") +
    xlab("") +
    theme_bw() +
    guides(fill = guide_legend(title = "Approach")) +
    # scale_fill_manual(values = ) +
    geom_col() +
    facet_wrap(~metric, ncol = 3, scales = "free") +
    geom_label(aes(label = round(value * 100, 2)), col = "white") +
    theme(axis.text.x = element_blank())
  pdf(fig_name, width = 6.5, height = 3.5)
  plot(p)
  dev.off()
  return(p)
}
plot_n_edges <- function(networks) {
  p_df <- data.frame(TYPE = names(networks), n_edges = sapply(networks, function(x) nrow(x$edges)))
  p <- ggplot(subset(p_df, TYPE != "GTRD"), aes(x = TYPE, y = n_edges, fill = TYPE)) +
    ylab("Number of edges") +
    xlab("") +
    theme_bw() +
    geom_col() +
    geom_label(aes(label = round(n_edges * 100, 2)), col = "white")
  plot(p)
}

get_tfs <- function(path_tf_and_reqdgenes) {
  input_tfs <- file.path(path_tf_and_reqdgenes, "tfs_anonymized.txt")
  read.table(input_tfs)[, 1]
}

read_graph_prediction <- function(TYPE, TEST_ID, cutoff) {
  # Load prediction results
  fname_prediction_scores <- paste0("../", TYPE, "/data/", TEST_ID, "_prediction_score.txt")
  df_prediction_scores <- read.table(fname_prediction_scores)
  colnames(df_prediction_scores) <- c("TF", "target", "score")
  df_prediction_scores_filt <- subset(df_prediction_scores, score > cutoff)
}
