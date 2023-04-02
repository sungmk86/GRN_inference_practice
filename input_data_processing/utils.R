library(igraph)
library(data.table)
library(ggplot2)
library(R6)
##########################

# Required functions for make_input.R
MakeInput <- R6Class("MakeInput",
    public = list(
        path_input = NULL,
        path_tf_and_reqdgenes = NULL,
        path_output = NULL,
        df_expr = NULL,
        df_net = NULL,
        tfs_all = NULL,
        required_genes = NULL,
        initialize = function(path_input, path_tf_and_reqdgenes, path_output, rng_seed = 1233) {
            set.seed(rng_seed)
            self$path_input <- path_input
            self$path_tf_and_reqdgenes <- path_tf_and_reqdgenes
            self$path_output <- path_output
            self$read_data()
        },
        read_data = function() {
            self$read_expression()
            self$read_tfs()
            self$read_required_genes()
            self$read_network()
        },
        read_metadata = function() {
            input_expr_meta <- file.path(
                self$path_input,
                "Bayesian_DE/iterative_test/iterative_test/TF_experiment_metadata.gz"
            )
            return(read.delim(input_expr_meta))
        },
        read_expression = function(cells_to_be_removed = c("iPSC", "K562"), column_name = "label.main", n_cells_for_pseudobulking = 150) {
            input_expr <- file.path(
                self$path_input,
                "Bayesian_DE/iterative_test/iterative_test/TF_experiment_expression_matrix.gz"
            )
            df_expr_full <- read.delim(input_expr, row.names = 1)
            df_expr_full$id <- NULL

            df_metadata <- self$read_metadata()
            message("Construct pseudobulk data")
            celltypes <- setdiff(df_metadata[[column_name]], cells_to_be_removed)
            list_pseudobulk <- sapply(celltypes, function(celltype) {
                message(celltype)
                selected_idx <- df_metadata[[column_name]] == celltype
                selected_ids <- rep(sample(colnames(df_expr_full)[selected_idx]),
                    length.out = n_cells_for_pseudobulking
                )
                apply(df_expr_full[, selected_ids], 1, sum)
            }, simplify = F)
            self$df_expr <- as.data.frame(list_pseudobulk)
        },
        read_tfs = function() {
            input_tfs <- file.path(self$path_tf_and_reqdgenes, "tfs_anonymized.txt")
            tfs <- read.table(input_tfs)[, 1]
            external_tfs <- grep("tf", rownames(self$df_expr), value = T)
            self$tfs_all <- unique(c(external_tfs, tfs))
        },
        read_required_genes = function() {
            input_required_genes <- file.path(self$path_tf_and_reqdgenes, "reqd_genes_anonymized.txt")
            self$required_genes <- read.table(input_required_genes)[, 1]
        },
        read_network = function() {
            input_network <- file.path(self$path_input, "BIC/data/networks_anonymize.txt")
            df_net <- read.delim(input_network)
            df_net <- df_net[!duplicated(df_net), ]
            colnames(df_net) <- c("TF", "target")
            self$df_net <- df_net
        },
        write_files = function(type) {
            selected_tfs <- intersect(rownames(self$df_expr), self$tfs_all)

            # write embeddings
            if (type == "predicting") {
                df_embeddings <- self$df_expr
                ###############################################
                # write all possible edges from TFs to target
                index_all_tfs <- match(selected_tfs, rownames(df_embeddings)) - 1
                index_all_genes <- seq_len(nrow(self$df_expr)) - 1
                all_possible_edges_index <- expand.grid(index_all_tfs, index_all_genes)
                write.table(all_possible_edges_index,
                    file.path(
                        self$path_output,
                        "edge_label_index_for_predicting.txt"
                    ),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )
                ###############################################
            } else if (type == "training") {
                selected_genes <- intersect(
                    rownames(self$df_expr),
                    c(self$tfs_all, self$required_genes)
                )
                df_embeddings <- self$df_expr[selected_genes, ]
            } else {
                stop("type can be either training or predicting.")
            }
            write.table(df_embeddings,
                paste0(self$path_output, "/embeddings_for_", type, ".txt"),
                quote = F, row.names = F, col.names = F, sep = "\t"
            )
            df_net_filt <- subset(
                self$df_net,
                TF %in% selected_tfs & target %in% rownames(df_embeddings)
            )
            df_edge_index <- data.frame(
                TF = match(df_net_filt$TF, rownames(df_embeddings)) - 1,
                target = match(df_net_filt$target, rownames(df_embeddings)) - 1
            )
            write.table(df_edge_index,
                paste0(self$path_output, "/edge_index_for_", type, ".txt"),
                quote = F, row.names = F, col.names = F, sep = "\t"
            )
        }
    )
)
