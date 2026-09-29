setwd(getwd())
library("tximport")
library("rhdf5")
library("DescTools")
library("DESeq2")
library("factoextra")

# Import the information about the sequenced samples
sinfo <- read.csv2("ref/sinfo.csv", header = TRUE, sep = ",")
profile <- data.frame(row.names = c(1:100))
for (i in list.files("scripts/batch_effect/qualimap/")) {
  fname <- paste0("scripts/batch_effect/qualimap/", i)
  sname <- substr(i,1,3)
  profile[sname] <- read.delim(fname)$Transcript.coverage.profile
}

# Normalize all the coverage values of each sample using the maximum value 
# from each sample
profile_maxnorm <- sweep(profile, 2, 
                         apply(profile, 2, max, na.rm = TRUE), FUN = "/")

# Create a trapezoid using the coverage values of each sample and calculate
# the surface values
x <- as.numeric(1:100)
sinfo["profile"] <- NA
for (i in 1:nrow(sinfo)) {
  sample <- as.character(sinfo[i, "samples"])
  print(sample)
  y <- profile_maxnorm[[sample]]
  sinfo[i, "profile"] <- AUC(x, y, method = "trapezoid")
}

# Clear the environment
rm(x, y, i, sample, profile, profile_maxnorm,fname,sname)

# Round the coverage surface values to 4 discrete groups
mx <- max(sinfo$profile)
mn <- min(sinfo$profile)
rn <- mx - mn
gr <- rn / 4
gr1 <- mx - gr
gr2 <- mx - (gr * 2)
gr3 <- mx - (gr * 3)
gr4 <- mx - (gr * 4)

sinfo["profile.r"] <- NA
for (i in 1:nrow(sinfo)) {
  if (sinfo[i, "profile"] > gr1) {
    sinfo[i, "profile.r"] <- 1
  } else if (sinfo[i, "profile"] > gr2) {
    sinfo[i, "profile.r"] <- 2
  } else if (sinfo[i, "profile"] > gr3) {
    sinfo[i, "profile.r"] <- 3
  } else if (sinfo[i, "profile"] > gr4) {
    sinfo[i, "profile.r"] <- 4
  } else {
    sinfo[i, "profile.r"] <- 4
  }
}
rm(mx, mn, rn, gr, gr1, gr2, gr3, gr4, i)
row.names(sinfo) <- sinfo$sample

# Import the transcript to gene file
tx <- read.table("ref/T2G.tsv", sep="\t", header = TRUE)

# Import the ribosomic gene list from gtf
ribo <- read.csv2("ref/ribo.txt", header = FALSE)

# Import kallisto quanification data
files <- file.path("kal_lnc", sinfo$sample, "abundance.h5")
names(files) <- paste0(sinfo$sample)
txi <- tximport(files, type = "kallisto", tx2gene = tx)

# Remove ribosomic genes from counts
counts <- as.data.frame(txi$counts)
counts <- counts[-which(row.names(counts) %in% ribo$V1), ]

#Clear the environment
rm(ribo, tx, files, files_out)

# Load sva package (Not before! DescTools and sva are not compatible)
library("sva")

# Correct counts with ComBat_seq
counts <- as.matrix(counts)
mode(counts) <- "numeric"
correctedCounts <- round(ComBat_seq(counts = counts, 
                                    batch = as.factor(sinfo$d_profile), 
                                    group = as.factor(sinfo$condition)))

# Import corrected data in DESeq2, normalize it and filter by expression
dds <- DESeqDataSetFromMatrix(countData = correctedCounts, 
                              colData = sinfo, design = ~ condition)
dds <- estimateSizeFactors(dds)
normCounts <- as.data.frame(counts(dds, normalized=TRUE))
fnormCounts <- normCounts[(rowSums(normCounts > 10)) > (length(normCounts)/2), ]
correctedCounts <- correctedCounts[rownames(correctedCounts) %in% 
                                     rownames(fnormCounts),]
dds <- DESeqDataSetFromMatrix(countData = correctedCounts, 
                              colData = sinfo, design = ~ condition)
dds <- estimateSizeFactors(dds)
normCounts <- counts(dds, normalized=TRUE)

# Perform Principal component analysis
pcaf <- prcomp(t(log(normcounts + 1)), scale. = FALSE, center = T)
p <- fviz_pca_ind(pcaf,
                  col.ind = sinfo$condition,
                  axes = c(1, 2),
                  addEllipses = FALSE,
                  ellipse.level = 0.95,
                  geom = c("point", "text"),
                  repel = TRUE
)
print(p)

# Perform the DESeq2 core function to identify DEGs
dds <- DESeq(dds)

# Select the p adjusted value and the log 2 fold change
pval <- 0.05
l2fc <- 0.58

# Extract Deseq2 results
dfresVC <- as.data.frame(results(dds_out,
  contrast=c("condition", "Vaccine", "Control"), alpha = pval))
dfresAC <- as.data.frame(results(dds_out,
  contrast=c("condition", "Adjuvant", "Control"), alpha = pval))
dfresVA <- as.data.frame(results(dds_out,
  contrast=c("condition", "Vaccine", "Adjuvant"), alpha = pval))

# Extract only significant DEGs (selected pval and l2fc)
sigVC <- as.data.frame(dfresVC[which(
  dfresVC$padj < pval & abs(dfresVC$log2FoldChange) > l2fc), ])
sigAC <- as.data.frame(dfresAC[which(
  dfresAC$padj < pval & abs(dfresAC$log2FoldChange) > l2fc), ])
sigVA <- as.data.frame(dfresVA[which(
  dfresVA$padj < pval & abs(dfresVA$log2FoldChange) > l2fc), ])

# Save data
save(dfresVC, dfresAC, dfresVA, sigVC, sigAC, sigVA, correctedCounts, 
     normCounts, dds, sinfo, file = "00_scripts/rdata/expression.RData")
