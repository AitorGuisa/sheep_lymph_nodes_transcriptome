setwd(getwd())
library("VennDiagram")
library("BioVenn")
library("eulerr")

# Load results from the three programs 
cpc2 <- read.table("lnc/CPC2_results_summary.txt", sep="\t")
cpat <- read.table("lnc/CPAT_results.txt", sep="\t")
hmmscan <- read.table("lnc/hmmscan_results.txt", sep="\t")

# Select only the non coding RNAs
ocpc2 <- cpc2[cpc2$V2 == "ncRNA", ]$V1
ocpat <- cpat[cpat$V2 == "ncRNA", ]$V1
ohmmscan <- hmmscan[hmmscan$V2 == "ncRNA", ]$V1

venn <- data.frame("intersections" = c("all", "only_cpat_cpc2", 
                                       "only_cpat_hmmer", "only_cpc2_hmmer", 
                                       "only_cpat", "only_cpc2", "only_hmmer"),
                   "values" = NA)

# Calculate the number of novel potenrial lncRNAs in each Venn diagram segment
venn[1,2] <- length(intersect(intersect(ocpat, ohmmscan), ocpc2))
venn[2,2] <- length(intersect(ocpat, ocpc2)) - venn[1,2]
venn[3,2] <- length(intersect(ocpat, ohmmscan)) - venn[1,2]
venn[4,2] <- length(intersect(ocpc2, ohmmscan)) - venn[1,2]
venn[5,2] <- length(ocpat) - venn[1,2] - venn[2,2] - venn[3,2]
venn[6,2] <- length(ocpc2) - venn[1,2] - venn[2,2] - venn[4,2]
venn[7,2] <- length(ohmmscan) - venn[1,2] - venn[3,2] - venn[4,2]

# plot venn diagram with eulerr (use values from venn)
fit <- euler(c("CPAT" = 281, "CPC2" = 234, "HMMER" = 929, "CPAT&CPC2" = 259, 
               "CPAT&HMMER" = 148, "CPC2&HMMER" = 1177, 
               "CPAT&CPC2&HMMER" = 1608))

lnc_venn <- plot(fit, quantities = TRUE)

# Write final lncRNA list and export it as "lncRNA_final_list.txt"
cpc2$Combined <- paste0(cpc2$V1,"-",cpc2$V2)
cpat$Combined <- paste0(cpat$V1,"-",cpat$V2)
hmmscan$Combined <- paste0(hmmscan$V1,"-",hmmscan$V2)
list <- intersect(cpc2$Combined,intersect(cpat$Combined,hmmscan$Combined))
data <- data.frame(matrix(unlist(strsplit(list,"-")), 
                          nrow=length(strsplit(list,"-")), byrow=TRUE))
sum(data$X2=="ncRNA")
write.table(data$X1[data$X2=="ncRNA"], file="lnc/lncRNA_final_list.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)
