setwd(getwd())
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(Rfuntzioak)
load("00_scripts/rdata/diff_expr.RData")

de_miRNAs <- unique(c(rownames(sigAC), rownames(sigVA), rownames(sigVC)))
rm(dfresAC, dfresVA, dfresVC, sigAC, sigVA, sigVC)

# Miranda: import data and filter by score and energy
miranda <- read.csv("05_target_prediction/ensembl/miranda/miranda.txt", 
                    sep = "\t")
miranda <- miranda %>% filter(score >= 145, energy < -20)
miranda_list <- list()
for (i in de_miRNAs) {
  miranda_list[[i]] <- unique(miranda[miranda$microRNA == i, "mRNA"])
}

# Pita: import data and filter by energy 
pita <- read.csv("05_target_prediction/ensembl/pita/pita.txt", sep = "\t")
pita <- pita %>% filter(energy < -20)
pita_list <- list()
for (i in de_miRNAs) {
  pita_list[[i]] <- unique(pita[pita$microRNA == i, "mRNA"])
}

# Risearch2: Import the data of each miRNA
risearch_list <- list()
for (i in de_miRNAs) {
  path <- "05_target_prediction/ensembl/risearch/"
  mir_file <- paste0(path, "risearch_", i, ".out.gz")
  risearch_list[[i]] <- as.vector(unique(
    read.table(file = mir_file, sep = "\t")$V4))
}

sum(lengths(miranda_list))
sum(lengths(pita_list))
sum(lengths(risearch_list))

# Intersect all the filtered target of each miRNA from the three tools
mi_ri_pi_intersect <- list()
for (i in de_miRNAs) {
  mi_ri_pi_intersect[[i]] <- intersect(intersect(miranda_list[[i]], 
                                        pita_list[[i]]), risearch_list[[i]])
}

# Clear enviroment
rm(miranda, pita, risearch_list, pita_list, miranda_list, i, mir_file, path)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTIONAL ANNOTATION
# Import the domain scope
load("00_scripts/rdata/ens_counts.RData")
domain_scope <- unique(rownames(normCounts))

mi_ri_pi_targets <- list()
for (i in de_miRNAs) {
  mi_ri_pi_targets[[i]] <-intersect(mi_ri_pi_intersect[[i]], domain_scope)
}

sum(lengths(mi_ri_pi_targets))

# This function converts sheep genes in human orthologs
ortol <- function(miRNA) {
  miRNA_ort <- gorth(
    miRNA,
    source_organism = "oarambouillet",
    target_organism = "hsapiens", 
    numeric_ns = "",
    mthreshold = Inf,
    filter_na = TRUE
  )
  miRNA_ort <- unique(miRNA_ort$ortholog_ensg)
  return(miRNA_ort)
}
domain_scope_hsa <- ortol(domain_scope)

# This function performs the ortol funcution and functional annotation analysis.
gostg_hsa <- function(miRNA, fname, inters) {
  mirna_ort <- ortol(miRNA)
  print(length(mirna_ort))
  gores <- gost(query = mirna_ort, 
                organism = "hsapiens", 
                evcodes = TRUE,
                correction_method = "fdr",
                user_threshold = 0.05,
                domain_scope = "custom", 
                custom_bg = domain_scope_hsa, 
                sources = NULL)
  gores <- apply(gores$result,2,as.character)
  gores <- gores[, colnames(gores) != "parents"]
  fpath <- paste0("07_func_annot/ensembl/", inters, "/hsa/", fname, ".csv")
  write.table(gores, fpath, sep = "\t", quote = FALSE, row.names = FALSE)
}

for (i in de_miRNAs) {
  print(i)
  try(gostg_hsa(mi_ri_pi_targets[[i]], i, "mi_ri_pi"))
}

save(mi_ri_pi_intersect, mi_ri_pi_targets,
     file = "00_scripts/rdata/miRNA_targets.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ENRICHMENT MAP ANALYSIS
# Perform functional annotation of all the DE miRNAs
all_genes <- c()
for (i in de_miRNAs) {
  mirgens <- mi_ri_pi_targets[[i]]
  all_genes <- c(all_genes, mirgens)
}
all_genes <- unique(all_genes)
all_genes_hsa <- ortol(all_genes)

gp <- gost(query = all_genes_hsa, 
           organism = "hsapiens",
           evcodes = TRUE,
           correction_method = "fdr",
           user_threshold = 0.05,
           domain_scope = "custom",
           significant = TRUE,
           custom_bg = domain_scope_hsa, 
           sources = c("GO:BP", "REAC", "KEGG"))
 
gp$result$intersection <- sapply(gp$result$intersection, 
                                 function(x) paste(x, collapse = ","))
extraxted_gp <- gp$result
extraxted_gp <- extraxted_gp[extraxted_gp["term_size"] > 5, ]
extraxted_gp <- extraxted_gp[extraxted_gp["term_size"] <= 500, ]
extraxted_gp <- extraxted_gp[,-14]
write.table(extraxted_gp, "07_func_annot/ensembl/mi_ri_pi/hsa/all_genes.csv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Create a table to input in Cytoscape to perform an EnrichmentMap analysis
res <- data.frame(
  GO.ID = extraxted_gp$term_id,
  Description = extraxted_gp$term_name,
  p.Val = extraxted_gp$p_value,
  FDR = extraxted_gp$p_value,
  Phenotype = "+1",
  Genes = extraxted_gp$intersection
)

# Write table
write.table(res, "08_cytoescape/all_gem.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
