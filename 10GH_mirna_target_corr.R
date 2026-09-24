setwd(getwd())
library("dplyr")
library("ggplot2")
library("gprofiler2")
load("00_scripts/rdata/miRNA_targets.RData")
load("00_scripts/rdata/expr.RData")
rnaseq_data <- "/media/labo/Expansion/lymph_node_total"
load(paste0(rnaseq_data, "/00_scripts/rdata/expression.RData"))

# Clear environment
rm(counts, correctedCounts, expressionMatrix, sinfo_1, 
   sinfo_out, dds_out, dfresAC, dfresVA, dfresVC)

# Rename data
miRNA_targets <- mi_ri_pi_intersect

# Filter target list to those that has three or more miRNA binding site
mult_bs_miranda <- 
  read.csv2("05_target_prediction/ensembl/miranda/multipleTargets_1_7.txt", 
            sep = "\t")
mult_bs_pita <- 
  read.csv2("05_target_prediction/ensembl/pita/multipleTargets_1_7.txt", 
            sep = "\t")
mult_bs_miranda <- mult_bs_miranda[mult_bs_miranda["numberTargets"] >= 3, ]
mult_bs_pita <- mult_bs_pita[mult_bs_pita["numberTargets"] >= 3, ]

# Intersect two data sets and split first column into mirnas and genes
mults_bs_intersect <- mult_bs_pita[mult_bs_pita$miRNA.transcript %in% 
                                     mult_bs_miranda$miRNA.transcript, ]
split_cols <- strsplit(mults_bs_intersect$miRNA.transcript, "#", fixed = TRUE)
mults_bs_intersect$mirna <- sapply(split_cols, `[`, 1)
mults_bs_intersect$gene <- sapply(split_cols, `[`, 2)

# Create a list 
miRNA_targets_mbs <- list()
for (i in unique(mults_bs_intersect$mirna)){
  vector <- mults_bs_intersect[mults_bs_intersect["mirna"] == i, ]
  miRNA_targets_mbs[[i]] <- as.vector(vector$gene)
}

# Select only that genes that has 3 or more binding sites and are that are 
# detected by tree all the programs
miRNA_targets_inters_mbs <- list()
for (i in names(miRNA_targets_mbs)) {
  genes_inter <- miRNA_targets[[i]]
  genes_mbs <- miRNA_targets_mbs[[i]]
  genes_inter_mbs <- genes_inter[genes_inter %in% genes_mbs]
  miRNA_targets_inters_mbs[[i]] <- genes_inter_mbs
}
miRNA_targets <- miRNA_targets_inters_mbs
sum(lengths(miRNA_targets))

# Convert sheep ensemble gene IDs to gene names
miRNA_targets_ncbi <- miRNA_targets
for (i in 1:length(miRNA_targets)) {
   print(i)
  x <- miRNA_targets[[i]]
  if (length(x) > 0){
    y <- gorth(x, "oarambouillet", "hsapiens")
    if (is.null(y)){
      miRNA_targets_ncbi[[i]] <- NA
    } else {
      miRNA_targets_ncbi[[i]] <- y$ortholog_name
    }
  } else {
    miRNA_targets_ncbi[[i]] <- NA
  }
}

# Obtain human orthologs (use the same function written in the script 05GH)
normCounts_out_noloc <- Sheep2Human(normCounts)

# Filter target list to those that are expressed in lymph nodes
miRNA_targets_expr <- miRNA_targets_ncbi
for (i in 1:length(miRNA_targets)) {
  targets <- miRNA_targets_ncbi[[i]]
  miRNA_targets_expr[[i]] <- 
    targets[targets %in% rownames(normCounts_out_noloc)]
}

# Create a data frame with all miRNA and target pairs
correlation_df <- data.frame()
for (i in 1:length(miRNA_targets_expr)) {
  print(i)
  miRNA <- names(miRNA_targets_expr)[i]
  miRNA <- gsub("_", "-", miRNA)
  target <- miRNA_targets_expr[[i]]
  tmp <- data.frame(
    miRNA = rep(miRNA, length(target)),
    target = target,
    stringsAsFactors = FALSE
  )
  correlation_df <- rbind(correlation_df, tmp)
}

# Keep only the samples that are present in both data sets
fnormCounts <- fnormCounts[ ,-c(6, 9)]
normCounts_out_noloc <- normCounts_out_noloc[ ,-6]

# Check if they match
colnames(fnormCounts)
colnames(normCounts_out_noloc)

# Perform correlation test between miRNA and target pairs
correlation_df["Correlation"] <- NA
correlation_df["p_value"] <- NA
for (i in 1:nrow(correlation_df)){
  mir <- correlation_df[i, "miRNA"]
  gen <- correlation_df[i, "target"]
  print(paste0(i, ": ", mir, " and ", gen))
  mir_e <- as.numeric(fnormCounts[mir, ])
  gen_e <- as.numeric(normCounts_out_noloc[gen, ])
  test <- cor.test(gen_e, mir_e, method = "spearman")
  correlation_df[i, "Correlation"] <- test$estimate
  correlation_df[i, "p_value"] <- test$p.value
}

# The pair novel-mir-47-m and SLC6A9 is duplicated, remove one!
correlation_df <- correlation_df[-28,]
correlation_df["p_adjusted"] <- p.adjust(correlation_df$p_value, method = "BH")
correlation_df <- correlation_df[order(correlation_df$p_adjusted), ]
write.table(correlation_df, "correlation.csv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Filter target list to those that are differentially expressed
hsigAC <- Sheep2Human(sigAC)
hsigVA <- Sheep2Human(sigVA)
hsigVC <- Sheep2Human(sigVC)
de_genes <- unique(c(rownames(hsigAC), rownames(hsigVA), rownames(hsigVC)))
for (i in 1:length(miRNA_targets)) {
  targets <- miRNA_targets_expr[[i]]
  miRNA_targets_expr[[i]] <- targets[targets %in% de_genes]
}
