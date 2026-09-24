setwd(getwd())
library("gprofiler2")
library("ggplot2")
library("rrvgo")

# Import the list of DEGs
load("00_scripts/rdata/expression.RData")

# Create a function to obtain human orthologs from sheep genes
Sheep2Human <- function(df) {
  df <- as.data.frame(df)
  loc_id <- df[startsWith(rownames(df), "LOC"), ]
  symbols <- df[!(rownames(df) %in% rownames(loc_id)), ]
  loc_ort <- gorth(
    rownames(loc_id),
    source_organism = "oarambouillet",
    target_organism = "hsapiens",
    numeric_ns = "",
    mthreshold = Inf,
    filter_na = TRUE)
  loc_id["symbol"] <- NA
  loc_id$symbol <- ifelse(rownames(loc_id) %in% loc_ort$input, 
    loc_ort$ortholog_name[match(rownames(loc_id), 
    loc_ort$input)], loc_id$symbol)
  loc_id <- loc_id[!is.na(loc_id$symbol), ]
  loc_id <- loc_id[!duplicated(loc_id$symbol), ]
  rownames(loc_id) <- loc_id$symbol
  loc_id <- subset(loc_id, select=-symbol)
  orthologs <- rbind(symbols, loc_id)
  return(orthologs)
}

# Obtain human orthologs for significant DEGs
hsigVC <- Sheep2Human(sigVC)
hsigAC <- Sheep2Human(sigAC)
hsigVA <- Sheep2Human(sigVA)
hdomain_scope <- Sheep2Human(normCounts_out)
hdomain_scope <- rownames(hdomain_scope)

# Remove the unannotated lncRNAs from domain_scope
lnc <- read.csv("00_scripts/classification/output_files/lncRNAs_gene.txt",
                sep = "\t", row.names = 1)
domain_scope <- rownames(normCounts_out)[! rownames(normCounts_out) 
                                         %in% rownames(lnc)] 

# Create a function to compute the functional annotation
gostg <- function(sig, org, file_name, dom) {
  gores <- gost(query = rownames(sig), 
                organism = org, 
                evcodes = TRUE,
                correction_method = "fdr",
                user_threshold = 0.05,
                domain_scope = "custom", 
                custom_bg = dom, 
                sources = NULL)
  gores <- apply(gores$result,2,as.character)
  gores <- gores[, colnames(gores) != "parents"]
  fpath <- paste0("out/func_annot/", file_name, ".csv")
  write.table(gores, fpath, sep = "\t", quote = FALSE, row.names = FALSE)
  return(gores)
}

hVC <- gostg(hsigVC, "hsapiens", "hVC", hdomain_scope)
hAC <- gostg(hsigAC, "hsapiens", "hAC", hdomain_scope)
hVA <- gostg(hsigVA, "hsapiens", "hVA", hdomain_scope)

# Perform the Treemap analisys
hva_df <- as.data.frame(hVA)
hvc_df <- as.data.frame(hVC)
hac_df <- as.data.frame(hAC)

bp_hva <- hva_df[hva_df["source"] == "GO:BP", ]
bp_hvc <- hvc_df[hvc_df["source"] == "GO:BP", ]
bp_hac <- hac_df[hac_df["source"] == "GO:BP", ]

simMatrix_hva <- calculateSimMatrix(bp_hva$term_id, 
                                    orgdb="org.Hs.eg.db", 
                                    ont="BP", 
                                    method="Rel")
simMatrix_hvc <- calculateSimMatrix(bp_hvc$term_id, 
                                    orgdb="org.Hs.eg.db", 
                                    ont="BP", 
                                    method="Rel")
simMatrix_hac <- calculateSimMatrix(bp_hac$term_id, 
                                    orgdb="org.Hs.eg.db", 
                                    ont="BP", 
                                    method="Rel")

scores_hva <- setNames(-log10(as.numeric(bp_hva$p_value)), bp_hva$term_id)
scores_hvc <- setNames(-log10(as.numeric(bp_hvc$p_value)), bp_hvc$term_id)
scores_hac <- setNames(-log10(as.numeric(bp_hac$p_value)), bp_hac$term_id)

reducedTerms_hva <- reduceSimMatrix(simMatrix_hva, 
                                    scores_hva,
                                    threshold=0.7, 
                                    orgdb="org.Hs.eg.db")                          
reducedTerms_hvc <- reduceSimMatrix(simMatrix_hvc, 
                                    scores_hvc,
                                    threshold=0.7, 
                                    orgdb="org.Hs.eg.db") 
reducedTerms_hac <- reduceSimMatrix(simMatrix_hac, 
                                    scores_hac,
                                    threshold=0.7, 
                                    orgdb="org.Hs.eg.db") 

# Save data
save(simMatrix_hva, simMatrix_hvc, reducedTerms_hva, reducedTerms_hvc,
     hVC, hAC, hVA, file = "00_scripts/rdata/func_annot.RData")
