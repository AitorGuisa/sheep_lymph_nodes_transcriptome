setwd(getwd())
library("ggplot2")
library("gprofiler2")
library("GWENA")
library("dplyr")
load("00_scripts/rdata/expression.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CO-EXPRESSION ANALYSIS
# Import the information about the sequenced samples.
sinfo["treatments"] <- c(rep("yestreated", 10) , rep( "notreated", 7))
normCounts_out <- as.data.frame(normCounts_out)
nc_filt <- filter_low_var(t(normCounts_out), pct = 0.7, type = "median")

# Build the networks. Set the R.sqrt threshold at 0.8.
net <- build_net(nc_filt, cor_func = "spearman", 
                 n_threads = 6, fit_cut_off = 0.8,
                 network_type = "signed")

# Extract the table of statistics for each power 
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

# Plot Scale independence and Mean connectivity scatter plots.
pow <- fit_power_table$Power
R2 <- fit_power_table$SFT.R.sq
mc <- fit_power_table$mean.k.
plot(pow, R2,
     type = "n",
     col = "blue",
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, Signed R²",
     main = "Scale Independence")
text(pow, R2, labels = pow, col = "blue", cex = 1)
abline(h = 0.8, col = "red", lty = 1, lwd = 1)
plot(pow, mc,
     type = "n",
     col = "blue",
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity")
text(pow, mc, labels = pow, col = "blue", cex = 1)

# Detect the modules in the networks
modules <- detect_modules(nc_filt, net$network, 
                          detailled_result = TRUE, 
                          merge_threshold = 0.9)

# GENES PER MODULE (BARPLOT)
# This commands create a barplot of the number of genes in each module.
genespermodule <- data.frame()
for (i in names(modules$modules)) {
  s <- data.frame(1)
  s["X1"] <- paste0("ME", i)
  s["genes"] <- length(modules$modules[[i]])
  genespermodule <- rbind(genespermodule, s)
} 
genespermodule$X1 <- factor(genespermodule$X1, 
                            levels = unique(genespermodule$X1))
ggplot(data=genespermodule, aes(x=X1, y=genes)) +
  geom_bar(stat="identity") +
  ylab("Number of genes") +
  xlab("Module")

# Calculate the correlation between modules and treatments
phenotype_association <- associate_phenotype(modules$modules_eigengenes, sinfo 
  %>% dplyr::select(condition, treatments))
plot_modules_phenotype(phenotype_association)
phenotype_association[["padj"]][["Vaccine"]] <- 
  p.adjust(phenotype_association$pval$Vaccine, method = "BH")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTIONAL ANNOTATION OF THE DETECTED MODULES
# Convert the sheep genes into human orthologs.
df <- as.data.frame(t(nc_filt))
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
  loc_ort$ortholog_name[match(rownames(loc_id), loc_ort$input)], loc_id$symbol)
loc_id <- loc_id[!is.na(loc_id$symbol), ]
loc_id <- loc_id[!duplicated(loc_id$symbol), ]
loc_id <- loc_id["symbol"]
symbols <- symbols["113"]
symbols["113"] <- rownames(symbols)
colnames(symbols) <- "symbol"
orthologs <- rbind(symbols, loc_id)

# Get the genes of each module and performs the functional annotation analysis.
nlist <- as.character(c(1:46))
module_func_prof <- list()
for (i in nlist) {
  print(i) 
  x <- orthologs[rownames(orthologs) %in% modules$modules[[i]], ]
  cat(as.character(length(x)), " orthologs in module ", i, "\n")
  gores <- tryCatch({
    gost(query = x, 
       organism = "hsapiens", 
       evcodes = TRUE,
       correction_method = "fdr",
       user_threshold = 0.05,
       domain_scope = "custom", 
       custom_bg = orthologs$symbol, 
       sources = NULL)
  }, error = function(e) {return(NULL)})
  
  if (is.null(gores) || is.null(gores$result) || nrow(gores$result) == 0) {
    message("no enrichment for module ", i)
    next
  }
  
  gores <- as.data.frame(apply(gores$result,2,as.character))
  gores <- gores[gores$source %in% c("GO:BP", "KEGG", "REAC"), ]
  title <- paste0("ME", as.character(i))
  module_func_prof[[title]] <- as.data.frame(gores)
}

func_annot_modules["query"] <- func_annot_modules$module
colnames(func_annot_modules)[1] <- "module"
func_annot_modules = func_annot_modules[func_annot_modules$source %in% 
                                          c("GO:BP", "KEGG", "REAC"), ]
write.table(func_annot_modules, "coexpr/func_annot_modules.csv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MODULE MEMBERSHIP
# This command return a data frame that shows the module membership of each gene
head(nc_filt)
aux <- as.data.frame(matrix(0,nrow=dim(nc_filt)[2],ncol=2))
colnames(aux) <- c("Gene","module")
aux$Gene <- colnames(nc_filt)
rownames(aux) <- aux$Gene
for(i in aux$Gene){
  for(j in names(modules$modules)){
    if(!is.na(match(i, modules$modules[[j]]))){
      aux[i,2] <- paste0("ME",j)
    }
  }
}
write.csv(aux, file="coexpr/module_membership.csv"
          ,quote=FALSE,col.names=TRUE,row.names=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LNCRNA STATISTICS
new_lnc <- read.delim("00_scripts/classification/output_files/lncRNAs_gene.txt",
                      row.names = 1)

de_new_vc <- sigVC[rownames(sigVC) %in% rownames(new_lnc), ]
de_new_ac <- sigAC[rownames(sigAC) %in% rownames(new_lnc), ]
de_new_va <- sigVA[rownames(sigVA) %in% rownames(new_lnc), ]

lnc_ms_vc <- aux[rownames(aux) %in% rownames(de_new_vc), ]
lnc_ms_ac <- aux[rownames(aux) %in% rownames(de_new_ac), ]
lnc_ms_va <- aux[rownames(aux) %in% rownames(de_new_va), ]

ms_vc <- aux[rownames(aux) %in% rownames(sigVC), ]
ms_ac <- aux[rownames(aux) %in% rownames(sigAC), ]
ms_va <- aux[rownames(aux) %in% rownames(sigVA), ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NODE AND EDGE TABLES FOR CYTOSCAPE
# Generate a node table with the gene amount per module and the association of 
# each module with the vaccine
node_table <- genespermodule
node_table["association_Vaccine"] <- phenotype_association$association$Vaccine
node_table["p-val_assoc_Vaccine"] <- phenotype_association$padj$Vaccine
write.table(node_table, "coexpr/node_tablen.csv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Generate a edge table
modules_eigengenes <- modules$modules_eigengenes
edge_table <- 
  as.data.frame(t(as.data.frame(combn(unique(colnames(modules_eigengenes)),2))))
colnames(edge_table) <- c("from", "to")
edge_table["R"] <- NA
edge_table["p-value"] <- NA
for (i in 1:nrow(edge_table)) {
  x <- edge_table[i, "from"]
  y <- edge_table[i, "to"]
  t <- cor.test(modules_eigengenes[[x]], modules_eigengenes[[y]])
  edge_table[i, "p-value"] <- t$p.value
  edge_table[i, "R"] <- as.character(t$estimate)
}
write.table(edge_table, "coexpr/edge_tablen.csv", sep = "\t", 
            quote = FALSE, row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Extract the hubgenes
eigenege <- modules$modules_eigengenes
gene_ModuleMembership <- as.data.frame(cor(nc_filt,eigenege, use="p"))
# Only for vaccine treatment 
gene_TraitSignificance <- as.data.frame(cor(
  nc_filt,c(rep(0,5), rep(1,5),rep(0,7)), use="p"))

# Select the vaccine related module, detect hub genes and perform enrichment
vac_rel_mod <- c("1","10")
hubgenes_func_prof <- list()
for (mod in vac_rel_mod) {
  # Extract the MM, GS and Names of the genes of the module (in order)
  MM <- gene_ModuleMembership[modules$modules[[mod]],as.numeric(mod)]
  GS <-gene_TraitSignificance[modules$modules[[mod]],]
  names <- modules$modules[[mod]]
  
  # Compute the 0.85 percentile of the MM and GS
  percentileMM <- quantile(MM,0.85)
  percentileGS <- quantile(abs(GS),0.85)
  
  # Select hub genes
  hubgenes <- names[MM>=percentileMM & abs(GS)>=percentileGS]
  
  # Perform functional annotation
  x <- orthologs[rownames(orthologs) %in% hubgenes, ]
  hub_enr <- gost(query = x, 
                  organism = "hsapiens", 
                  evcodes = TRUE,
                  correction_method = "fdr",
                  user_threshold = 0.05,
                  domain_scope = "custom", 
                  custom_bg = orthologs$symbol, 
                  sources = NULL)
  hub_enr <- as.data.frame(apply(hub_enr$result,2,as.character))
  hub_enr <- hub_enr[hub_enr$source %in% c("GO:BP", "KEGG", "REAC"), ]
  # Create a nested list
  title <- paste0("ME", as.character(mod))
  hubgenes_func_prof[[title]] <- list(hub_genes = hubgenes,
                                      fun_prof = hub_enr)

}
save(hubgenes_func_prof, module_func_prof, aux, 
     file = "scripts/rdata/hub_mod_enrich.RData")
