source("plot_utils.R")

tfnames <- read.table("../network/data/human_TFs_v2.7_mapped_to_Ensembl98.txt")[, 1]

df_deg <- read.table("../fantom5/data/deg_dermal.fibroblast_iPSC.txt")
norm_factor <- max(abs(df_deg$stat))
df_deg$stat_norm <- df_deg$stat / norm_factor

get_centrality <- function(df_net_weight, type = "default") {
    if (type == "default") {
        df_net <- df_net_weight[, c("from", "to")]
        g <- igraph::graph_from_data_frame(df_net, directed = T)
    } else if (type == "weight") {
        g <- igraph::graph_from_data_frame(df_net_weight, directed = T)
    } else if (type == "subset") {
        cutoff <- quantile(df_net_weight$weight, .3)
        g <- igraph::graph_from_data_frame(subset(df_net_weight, weight > cutoff)[,1:2], directed = T)
    }
    centrality <- eigen_centrality(g)$vector
    return(centrality)
}

lst_res <- lst_score <- list()
# for (TESTID in c("TEST15", "TEST15v2", "TEST15v3", "TEST15v4", "TEST15v5")) {
for (TESTID in c("TEST15v3", "TEST15v4", "TEST15v5")) {
    print(TESTID)
    df_combined_edge_scores <- combine_rng_edge_scores(TESTID, 50)
    lst_res[[TESTID]] <- apply(df_combined_edge_scores, 1, median)

    df_net <- matrix(unlist(strsplit(names(lst_res[[TESTID]]), "_")), byrow = T, ncol = 2)
    g_weight <- igraph::graph_from_data_frame(df_net, directed = T)
    E(g_weight)$weight <- lst_res[[TESTID]]
    df_net_weight <- as_data_frame(g_weight, what = "edges")

    df_out <- ddply(df_net_weight, "from", summarise, n_target = length(weight), weight_sum = sum(weight))
    df_out$cent <- get_centrality(df_net_weight, "default")[df_out$from]
    df_out$weightcent <- get_centrality(df_net_weight, "weight")[df_out$from]
    df_out$subsetcent <- get_centrality(df_net_weight, "subset")[df_out$from]

    df_out_deg <- merge(df_out, df_deg, by.x = "from", by.y = 0)
    df_out_deg_subset <- merge(df_out_subset, df_deg, by.x = "from", by.y = 0)
    df_out_deg <- transform(df_out_deg,
        score_cent = cent,
        score_weightcent = weightcent,
        score_subsetcent = subsetcent,
        score_cent_stat = cent * stat_norm,
        score_weightcent_stat = weightcent * stat_norm,
        score_subsetcent_stat = subsetcent * stat_norm
    )
    lst_score[[TESTID]] <- subset(df_out_deg, from %in% tfnames)
}

mat <- matrix("-", 20, 19)
for (idx in 1:19) {
    if (idx %in% 1:7) {
        df <- lst_score[["TEST15v5"]]
    } else if (idx %in% 8:13) {
        df <- lst_score[["TEST15v3"]]
    } else if (idx %in% 14:19) {
        df <- lst_score[["TEST15v4"]]
    }
    if (idx == 1) {
        df <- df[order(df$stat_norm, decreasing = T), ]
    } else if (idx %in% c(2, 8, 14)) {
        df <- df[order(df$score_cent, decreasing = T), ]
    } else if (idx %in% c(3, 9, 15)) {
        df <- df[order(df$score_weightcent, decreasing = T), ]
    } else if (idx %in% c(4, 10, 16)) {
        df <- df[order(df$score_subsetcent, decreasing = T), ]
    } else if (idx %in% c(5, 11, 17)) {
        df <- df[order(df$score_cent_stat, decreasing = T), ]
    } else if (idx %in% c(6, 12, 18)) {
        df <- df[order(df$score_weightcent_stat, decreasing = T), ]
    } else if (idx %in% c(7, 13, 19)) {
        df <- df[order(df$score_subsetcent_stat, decreasing = T), ]
    }
    mat[, idx] <- head(df[, c("from"), drop = F], 20)$from
}
