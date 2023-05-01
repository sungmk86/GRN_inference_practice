library(R6)
##########################

# Required functions for make_input.R
MakeInput <- R6Class("MakeInput",
    public = list(
        TEST_ID = NULL,
        path_tf_and_reqdgenes = NULL,
        path_output = NULL,
        df_expr = NULL,
        df_net = NULL,
        tfs_all = NULL,
        # required_genes = NULL,
        initialize = function(TEST_ID = NULL,
                              path_expr = NULL, path_meta = NULL,
                              pseudobulking = T, n_cells_for_selecting = 150, column_name = "label.main",
                              is_normalized = F, is_scaled = F,
                              path_network = NULL, path_tf_and_reqdgenes = NULL,
                              path_output = NULL, rng_seed = 1234) {
            set.seed(rng_seed)
            self$TEST_ID <- TEST_ID
            self$path_tf_and_reqdgenes <- path_tf_and_reqdgenes
            self$path_output <- path_output

            # read data
            self$read_expression(path_expr, path_meta, pseudobulking, n_cells_for_selecting, column_name, is_normalized, is_scaled)
            self$read_tfs()
            # self$read_required_genes()
            self$read_network(path_network)
        },
        read_expression = function(path_expr, path_meta, pseudobulking, n_cells_for_selecting, column_name, is_normalized, is_scaled, cells_to_be_removed = c("iPSC", "K562")) {
            df_expr_full <- read.delim(path_expr, row.names = 1)
            df_expr_full$id <- NULL

            df_metadata <- read.delim(path_meta)
            celltypes <- setdiff(df_metadata[[column_name]], cells_to_be_removed)

            if (pseudobulking) {
                message("Construct pseudobulk data")
                list_selected <- sapply(celltypes, function(celltype) {
                    message(celltype)
                    selected_idx <- df_metadata[[column_name]] == celltype
                    selected_ids <- rep(sample(colnames(df_expr_full)[selected_idx]),
                        length.out = n_cells_for_selecting
                    )
                    apply(df_expr_full[, selected_ids], 1, sum)
                }, simplify = F)
            } else {
                list_selected <- sapply(celltypes, function(celltype) {
                    message(celltype)
                    selected_idx <- df_metadata[[column_name]] == celltype
                    selected_ids <- rep(sample(colnames(df_expr_full)[selected_idx]),
                        length.out = n_cells_for_selecting
                    )
                    df_expr_full[, selected_ids]
                }, simplify = F)
            }
            df_combined <- as.data.frame(list_selected)
            # Filter genes
            df_combined_subset <- df_combined[apply(df_combined != 0, 1, sum) != 0, ]
            # Normalize
            if (!is_normalized) {
                df_norm <- t(t(df_combined_subset) / apply(df_combined_subset, 2, sum))
            } else {
                df_norm <- df_combined_subset
            }
            # Scaling
            if (!is_scaled) {
                df_scaled <- t(apply(log2(df_norm + 1), 1, scale))
            } else {
                df_scaled <- df_norm
            }
            colnames(df_scaled) <- colnames(df_combined_subset)
            self$df_expr <- df_scaled
        },
        read_tfs = function() {
            if (file.exists(self$path_tf_and_reqdgenes)) {
                input_tfs <- "/home/seongwonhwang/Desktop/projects/git/Node_ablation_practice/TFDB3/Gallus_gallus_TF"
                tfs <- setdiff(read.delim(input_tfs)$Symbol, "-")
                self$tfs_all <- tfs
            } else {
                input_tfs <- file.path(self$path_tf_and_reqdgenes, "tfs_anonymized.txt")
                tfs <- read.table(input_tfs)[, 1]
                external_tfs <- grep("tf", rownames(self$df_expr), value = T)
                self$tfs_all <- unique(c(external_tfs, tfs))
            }
        },
        read_required_genes = function() {
            input_required_genes <- file.path(self$path_tf_and_reqdgenes, "reqd_genes_anonymized.txt")
            self$required_genes <- read.table(input_required_genes)[, 1]
        },
        read_network = function(path_network) {
            df_net <- read.delim(path_network)
            df_net <- df_net[!duplicated(df_net), ]
            colnames(df_net) <- c("TF", "target")
            self$df_net <- df_net
        },
        write_files = function(type) {
            # selected_tfs <- intersect(rownames(self$df_expr), self$tfs_all)

            genes_in_network <- unique(c(self$df_net[, 1], self$df_net[, 2]))
            candidates <- sort(intersect(genes_in_network, rownames(self$df_expr)))

            # selected_genes <- intersect(
            #     rownames(self$df_expr),
            #     c(self$tfs_all, self$required_genes)
            # )
            # selected_genes <- rownames(self$df_expr)
            selected_genes <- candidates
            df_embeddings <- self$df_expr[selected_genes, ]
            g_cand <- subset(self$df_net, TF %in% intersect(candidates, self$tfs_all) & target %in% candidates)
            # write embeddings
            if (type == "predicting") {
                ###############################################
                # write all possible edges from TFs to target
                # index_all_tfs <- match(selected_tfs, rownames(df_embeddings)) - 1
                # index_all_genes <- seq_len(nrow(df_embeddings)) - 1
                # df_all_possible_edges_idx <- expand.grid(index_all_tfs, index_all_genes)

                df_net_filt <- g_cand
                df_all_possible_edges_idx <- data.frame(
                    TF = match(df_net_filt$TF, rownames(df_embeddings)) - 1,
                    target = match(df_net_filt$target, rownames(df_embeddings)) - 1
                )

                write.table(df_all_possible_edges_idx,
                    paste0(
                        self$path_output, "/",
                        self$TEST_ID, "_edge_index_for_predicting.txt"
                    ),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )
                # df_all_possible_edges <- expand.grid(selected_tfs, rownames(df_embeddings))

                df_all_possible_edges <- g_cand
                write.table(df_all_possible_edges,
                    paste0(
                        self$path_output, "/",
                        self$TEST_ID, "_all_possible_edges.txt"
                    ),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )
                ###############################################
            } else if (type == "training") {
                # df_net_filt <- subset(
                #     self$df_net,
                #     TF %in% selected_tfs & target %in% rownames(df_embeddings)
                # )

                df_net_filt <- g_cand
                df_edge_index <- data.frame(
                    TF = match(df_net_filt$TF, rownames(df_embeddings)) - 1,
                    target = match(df_net_filt$target, rownames(df_embeddings)) - 1
                )
                write.table(df_edge_index,
                    paste0(self$path_output, "/", self$TEST_ID, "_edge_index_for_", type, ".txt"),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )
            }
            write.table(df_embeddings,
                paste0(self$path_output, "/", self$TEST_ID, "_embeddings_for_", type, ".txt"),
                quote = F, row.names = F, col.names = F, sep = "\t"
            )
        }
    )
)
