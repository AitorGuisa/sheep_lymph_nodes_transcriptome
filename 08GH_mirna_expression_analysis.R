setwd(getwd())
library("dplyr")
library("DESeq2")
library("edgeR")
library("seqinr")
library("factoextra")

# Import quantification data and sample information data frames
counts <- read.csv('04_quantification/miRNAs_expressed_all_samples_now.csv',
                   sep="\t")
sinfo <- read.csv('00_references/sample_info.csv')
rownames(sinfo) <- sinfo$name

# Remove unnecessary columns from sample info and quantification data frames
sinfo <- sinfo[,-c(1,2)]
counts <- select(counts, -read_count, -precursor, -total)
counts <- select(counts, -X001.norm., -X002.norm., -X003.norm., 
                       -X004.norm., -X005.norm., -X006.norm., -X007.norm., 
                       -X008.norm., -X009.norm., -X010.norm., -X011.norm., 
                       -X012.norm., -X013.norm., -X014.norm., -X015.norm., 
                       -X016.norm., -X017.norm., -X018.norm., -X019.norm.)
 
# Remove duplicated rows
dup <- counts[duplicated(counts$X.miRNA),]
counts <- counts[-as.numeric(rownames(dup)), ]
rm(dup)
sinfo["sname"] <- paste0("ln_", sinfo$name)

# Set mirnas as rownames and samples as colnames
rownames(counts) <- counts[,1]
counts <- counts[,-1]
colnames(counts) <- sinfo$sname
colnames(counts) <- sinfo$name

# Normalize the data with CPM and filter by expression
normCounts <- as.data.frame(cpm(counts_1))
fnormCounts <- normCounts[(rowSums(normCounts > 1)) > (length(normCounts)/2), ]
expressionMatrix <- counts_1[rownames(counts_1) %in% rownames(fnormCounts), ]

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

# Remove potential outlier samples
sinfo <- sinfo[-6, ]
counts <- counts[,-6]

# Set adjusted p value and log2 fold change
pval <- 0.05
l2fc <- 0.58

# Perform differential expression analysis with DESeq2 and extract results
dds <- DESeqDataSetFromMatrix(countData = expressionMatrix, 
                              colData = sinfo, 
                              design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

dfresVC <- as.data.frame(results(dds, 
  contrast=c("condition","Vaccine", "Control"), alpha = pval))
dfresAC <- as.data.frame(results(dds,
  contrast=c("condition","Adjuvant", "Control"), alpha = pval))
dfresVA <- as.data.frame(results(dds,
  contrast=c("condition","Vaccine", "Adjuvant"), alpha = pval))

sigVC <- as.data.frame(dfresVC[which(dfresVC$padj < pval & 
  abs(dfresVC$log2FoldChange) > l2fc), ])
sigAC <- as.data.frame(dfresAC[which(dfresAC$padj < pval & 
  abs(dfresAC$log2FoldChange) > l2fc), ])
sigVA <- as.data.frame(dfresVA[which(dfresVA$padj < pval & 
  abs(dfresVA$log2FoldChange) > l2fc), ])

# save data
save(dfresAC, dfresVA, dfresVC, sigAC, sigVA, sigVC, 
     file = "00_scripts/rdata/diff_expr1.RData")

# Import all the miRNA sequences and write the significantly differentially 
# expressed mirna sequences in fasta files
miRNA_sequences <- read.fasta("00_references/total_seq.fa")
selected_mirnas <- unique(c(rownames(sigAC), rownames(sigVA), rownames(sigVC)))
miRNA_sequences <- miRNA_sequences[names(miRNA_sequences) %in% selected_mirnas]
write.fasta(sequences = miRNA_sequences, names = names(miRNA_sequences), 
            file.out = "05_target_prediction/miRNA_sequences.fa")







