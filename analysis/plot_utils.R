library(igraph)
library(ggplot2)
library(data.table)
library(cowplot)
library(ggrepel)
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
  path_networks <- c(
    list_predicted,
    list(
      eBIC = "/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/BIC/edge_strength_filt.txt",
      STRING = "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/BIC/data/networks_anonymize.txt",
      GTRD = "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/gtrd_anonymized.txt",
      Dorothea = "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/dorothea_anonymized.txt"
    )
  )

  networks <- sapply(path_networks, function(x) read_network(x, tfs_all, directed = T), simplify = F)
  return(networks)
}


plot_f1 <- function(results_f1, TEST_ID, type) {
  fig_name <- paste0("data/", TEST_ID, "_", type, ".pdf")
  accuracies <- sapply(results_f1, function(x) x$byClass)
  y <- melt(data.frame(metric = rownames(accuracies), accuracies))

  # p <- ggplot(
  #   y[y$metric %in% c("F1", "Recall", "Precision") & y$variable != type, ],
  #   aes(x = variable, y = value * 100, fill = variable)
  # ) +
  #   ylab("(%)") +
  #   xlab("") +
  #   theme_bw() +
  #   guides(fill = guide_legend(title = "Approach")) +
  #   # scale_fill_manual(values = ) +
  #   geom_col() +
  #   facet_wrap(~metric, ncol = 3, scales = "free") +
  #   geom_label(aes(label = round(value * 100, 2)), col = "white") +
  #   theme(axis.text.x = element_blank())
  # fig_width <- max(6.5, length(results_f1) * 1.5)
  p_df = y[y$metric %in% c("F1", "Recall", "Precision") & y$variable != type, ]
  p_df$group = sub('_rng[0-9]*', '', p_df$variable)
  p <- ggplot(
    p_df,
    aes(x = group, y = value * 100)
  ) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
    facet_wrap(~metric, ncol = 3, scales = "free") +
    ylab("(%)") +
    xlab("")
  fig_width <- max(4, length(results_f1) *.6)
  pdf(fig_name, width = fig_width, height = 3.5)
  plot(p)
  dev.off()
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

read_graph_prediction <- function(TYPE, TEST_ID, cutoff = NULL) {
  # Load prediction results
  fname_prediction_scores <- paste0("../", TYPE, "/data/", TEST_ID, "_prediction_score.txt")
  df_prediction_scores <- read.table(fname_prediction_scores)
  colnames(df_prediction_scores) <- c("TF", "target", "score")
  if (is.null(cutoff)) cutoff <- quantile(df_prediction_scores$score, .7)
  df_prediction_scores_filt <- subset(df_prediction_scores, score > cutoff)
}



plot_score_distribution <- function(dir_network, TEST_ID) {
  lst_predictions <- list.files(path = dir_network, pattern = "_prediction_score.txt")
  lst_predictions <- grep(TEST_ID, lst_predictions, value = T)

  lst_plots <- list()
  for (fname_pred in lst_predictions) {
    ID <- sub("_prediction_score.txt", "", fname_pred)
    df_pred <- read.table(file.path(dir_network, fname_pred), sep = "\t")

    colnames(df_pred) <- c("source", "target", "score")
    cutoff <- quantile(df_pred$score, .7)
    p <- ggplot(df_pred, aes(score)) +
      geom_histogram() +
      ggtitle(ID, subtitle = paste0("cutoff: ", round(cutoff, 2))) +
      geom_vline(xintercept = cutoff, lty = 2, colour = "red") +
      theme_classic() +
      coord_cartesian(xlim = c(0, 1))
    lst_plots[[ID]] <- p
  }
  fig_height <- max(length(lst_plots) / 3 * 500, 1000)
  png(paste0("data/score_distribution_", TEST_ID, ".png"), width = 1400, height = fig_height, res = 150)
  p <- plot_grid(plotlist = lst_plots, ncol = 2)
  plot(p)
  dev.off()
}


plot_scores <- function(df_scores, xlabel, ylabel, title, show_genenames = T) {
  p <- ggplot(df_scores, aes(x = score_1, y = score_2))

  if ("is_TF" %in% colnames(df_scores)) {
    p <- p + geom_point(aes(colour = is_TF), size = 2, alpha = .6) +
      scale_colour_manual(values = c("grey10", "red"))
  } else {
    p <- p + geom_point(size = 2, alpha = .6)
  }
  if (show_genenames) {
    # geom_text_repel(
    #   data = subset(df_scores, gene %in% tfnames & score_diff >= sort(subset(df_scores, gene %in% tfnames)$score_diff, decreasing = T)[3]),
    #   aes(label = gene), size = 4, nudge_y = .15
    # ) +
    p <- p +
      geom_text_repel(
        data = subset(df_scores, score_diff >= sort(df_scores$score_diff, decreasing = T)[3]),
        aes(label = gene), size = 4, nudge_y = -.15
      ) +
      geom_text_repel(
        data = subset(df_scores, score_diff <= sort(df_scores$score_diff, decreasing = F)[3]),
        aes(label = gene), size = 4, nudge_y = .15
      )
  }
  p <- p +
    theme_classic() +
    geom_abline(slope = 1, intercept = 0, lty = 2, colour = "red") +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title)

  return(p)
}

build_igraph <- function(fname_network) {
  df_net <- read.table(fname_network, sep = "\t")
  colnames(df_net) <- c("TF", "target", "score")

  # df_net_subset <- subset(df_net, TF %in% tfnames, select = c("TF", "target"))
  g <- igraph::graph_from_data_frame(df_net, directed = T)
  E(g)$weight <- df_net$score

  return(list(g = g, net = df_net))
}


compare_scores <- function(g1, g2, tfnames, type) {
  if (type == "eigen_centrality") {
    lst_scores1 <- eigen_centrality(g1, weights = E(g1)$weights)$vector
    lst_scores2 <- eigen_centrality(g2, weights = E(g2)$weights)$vector
  } else if (type == "degree") {
    lst_scores1 <- degree(g1)
    lst_scores2 <- degree(g2)
  }

  df_scores <- merge(
    data.frame(gene = names(lst_scores1), score = lst_scores1),
    data.frame(gene = names(lst_scores2), score = lst_scores2),
    by = "gene", suffixes = c("_1", "_2")
  )

  df_scores <- transform(df_scores,
    is_TF = gene %in% tfnames,
    score_diff = score_1 - score_2
  )
  return(df_scores)
}

plot_scores_comp <- function(
    dir_network, TEST_ID, lst_rng_seed, type, tfnames) {
  lst_plots <- list()
  for (idx1 in 1:(length(lst_rng_seed) - 1)) {
    for (idx2 in idx1:length(lst_rng_seed)) {
      rng_seed1 <- lst_rng_seed[idx1]
      rng_seed2 <- lst_rng_seed[idx2]

      if (idx1 == idx2) next
      ID <- paste0(sub("^_", "", rng_seed1), "_vs_", sub("^_", "", rng_seed2))
      print(ID)
      net1 <- build_igraph(paste0(dir_network, "/", TEST_ID, rng_seed1, "_prediction_score.txt"))
      net2 <- build_igraph(paste0(dir_network, "/", TEST_ID, rng_seed2, "_prediction_score.txt"))

      df_scores <- compare_scores(net1$g, net2$g, tfnames, type)
      lst_plots[[ID]] <- plot_scores(df_scores,
        xlabel = sub("^_", "", rng_seed1), ylabel = sub("^_", "", rng_seed2),
        title = paste0(type, "\n", ID)
      )
    }
  }
  fig_height <- max(length(lst_plots) / 3 * 500, 1300)
  png(paste0("data/", type, "_comp_between_rng_", TEST_ID, ".png"), width = 2100, height = fig_height, res = 150)
  p <- plot_grid(plotlist = lst_plots, ncol = 3)
  plot(p)
  dev.off()
}

plot_edge_scores_comp <- function(
    dir_network, TEST_ID, lst_rng_seed) {
  lst_plots <- list()
  for (idx1 in 1:(length(lst_rng_seed) - 1)) {
    for (idx2 in idx1:length(lst_rng_seed)) {
      rng_seed1 <- lst_rng_seed[idx1]
      rng_seed2 <- lst_rng_seed[idx2]

      if (idx1 == idx2) next
      ID <- paste0(sub("^_", "", rng_seed1), "_vs_", sub("^_", "", rng_seed2))
      print(ID)
      net1 <- build_igraph(paste0(dir_network, "/", TEST_ID, rng_seed1, "_prediction_score.txt"))
      net2 <- build_igraph(paste0(dir_network, "/", TEST_ID, rng_seed2, "_prediction_score.txt"))

      df_scores <- merge(net1$net, net2$net, by = c("TF", "target"), suffixes = c("_1", "_2"))
      cor_val <- with(df_scores, cor(score_1, score_2))
      lst_plots[[ID]] <- plot_scores(df_scores,
        xlabel = sub("^_", "", rng_seed1), ylabel = sub("^_", "", rng_seed2),
        title = paste0("Edge scores\n", ID, "\ncor: ", round(cor_val, 2)), show_genenames = F
      )
    }
  }
  fig_height <- max(length(lst_plots) / 3 * 500, 1300)
  png(paste0("data/Edge_scores_comp_between_rng_", TEST_ID, ".png"), width = 2100, height = fig_height, res = 150)
  p <- plot_grid(plotlist = lst_plots, ncol = 3)
  plot(p)
  dev.off()
}
