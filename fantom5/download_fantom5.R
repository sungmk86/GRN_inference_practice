# wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v7/extra/CAGE_peaks_df_fantom5ession//hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz -O data/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz
# gunzip data/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz
fname <- "data/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt"
df_fantom5 <- read.table(fname, sep = "\t", check.names = F, row.names = 1, header = T)
df_fantom5 <- df_fantom5[!grepl("STAT:", rownames(df_fantom5)), ]
df_fantom5 <- subset(df_fantom5, !is.na(short_description) & grepl("@", short_description) & !grepl("chr[0-9MX]*:|chrY:", short_description))

info <- strsplit(df_fantom5$short_description, ",")
df_fantom5 <- df_fantom5[rep(1:nrow(df_fantom5), sapply(info, length)), ]

geneSymbol <- toupper(gsub("^p[0-9]*@", "", unlist(info)))
df_expr <- apply(df_fantom5[, -c(1:6)], 2, tapply, geneSymbol, sum)

names_expr <- sapply(colnames(df_expr), function(x) sub(".hg38.nobarcode$", "", URLdecode(x)))
ffid <- sapply(strsplit(names_expr, "\\."), function(x) x[length(x)])
riken_id <- sapply(strsplit(names_expr, "\\."), function(x) tolower(x[length(x) - 1]))
description <- sapply(strsplit(names_expr, "\\."), function(x) paste(x[2:(length(x) - 2)], collapse = "."))

colnames(df_expr) <- riken_id
df_meta <- data.frame(ffid, riken_id, description)
write.table(df_expr, "data/hg38_fantom_tpm.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(df_meta, "data/hg38_fantom_meta.txt", quote = F, row.names = F, col.names = T, sep = "\t")


df_conv_sample_ids = read.table('data/hg38_conv_sample_ids.txt')
df_conv_list = read.table('data/hg38_conv_list.txt', fill =T)

conv_id = 'fib_to_H9ES'
