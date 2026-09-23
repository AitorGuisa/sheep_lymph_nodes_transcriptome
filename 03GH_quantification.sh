#!/bin/bash

base="/media/labo/Expansion/lymph_node_total"

# Using the list of concordant lncRNAs (lncRNA_final_list.txt), classify and name them according to position to nearly 
# annotated protein coding genes. We use bedtools (v2.25.0) for some steps, specifically bedtools closest module.
# First, we create a gtf file with the lncRNAs that passed al the filtering steps
cd $base/lnc
awk -F "\t" 'BEGIN{OFS="\t"} FNR==NR {myarray[$1]++; next} {aux=$9; gsub(/^transcript_id \"/,"",aux); gsub(/\".+/,"",aux); if( myarray[aux] ) {print $0}}' \
	lncRNA_final_list.txt lncRNA_candidate.gtf > lncRNA_final_candidate.gtf

# Check if some transcript classified as intronic (i) by gffcompare, its located really in other strands.
# If so, change its category to antisense (x).
# First create a file with the strand information of annotated transcripts
annotation="/media/labo/Expansion/lymph_node_total/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.gtf"
awk -F "\t" 'BEGIN{OFS="\t"} $3=="transcript" {gsub(/.*transcript_id \"/,"",$9); gsub(/\";.+/,"",$9); \
print $9,$7}' "$annotation" > transcript_strand_info.txt

# Change classification of intronic transcript located in the other strand
awk -F "\t" 'BEGIN{OFS="\t"} FNR==NR {strand[$1]=$2; next} $9 !~ /class_code \"i\"/ {print $0} \
$9 ~ /class_code \"i\"/ {gene=$9; gsub(/.*cmp_ref \"/,"",gene); gsub(/\";.+/,"",gene); \
if(strand[gene]==$7) {print $0} else {gsub(/\"i\"/,"\"x\"",$0); print $0}}' \
transcript_strand_info.txt lncRNA_final_candidate.gtf \
> lncRNA_final_candidate2.gtf && mv -f lncRNA_final_candidate2.gtf lncRNA_final_candidate.gtf

# Apply classification and naming algorithm (in temporary folder at home)
./classification.sh --gtf1 lncRNA_final_candidate.gtf --gtf2 $base/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic_nogenes.gtf -d 5000

# KALLISTO: We generate the indexes for the quantification and we introduce the index and the trimmed samples in kallisto to quantify the transcripts. 
#Sequence strandness was previously checked, to confirm that sequences are strand-specific and reverse strand have to be read firstly. 
#The indexex were constructed using the ARS-UI_Ramb_v2.0 and ARS-UI_Ramb_v3.0 sheep transcriptomes.
kallisto index --index=ref/genome_index/ind_v3_and_lnc.idx 00_scripts/classification/output_files/annotationANDlncRNAs.fa \
	2> ref/genome_index/ind_v3_and_lnc.log
cd rless
for f1 in *_1_p.fastq.gz
do
     f2=${f1%%_1_p.fastq.gz}"_2_p.fastq.gz"
     f3=${f1%%_1_p.fastq.gz}
     f4=${f1%%_1_p.fastq.gz}".log"
echo "quantifying $f1 and $f2 in $f3"
cd $base
mkdir -p kal_lnc/$f3
kallisto quant --rf-stranded --bias -t 8 -i ref/genome_index/ind_v3_and_lnc.idx -o kal_lnc/$f3 rless/$f1 rless/$f2 2> kal_lnc/$f3/$f4
done

# Rename directory names
cd $base/kal_lnc
while IFS=',' read -r original_name new_name; do
    mv "$original_name" "$new_name"
done < $base/00_scripts/directory_names.csv

# For quantification at gene level, create a transcript to gene file
echo -e "TXNAME\tGENEID" > T2G.tsv
awk -F "\t" 'BEGIN{OFS="\t"} $3=="transcript" {t_ID=$9; gsub(/.*transcript_id \"/,"",t_ID); gsub(/\";.+/,"",t_ID); \
g_ID=$9; gsub(/.*gene_id \"/,"",g_ID); gsub(/\";.+/,"",g_ID); print t_ID,g_ID \
}' 00_scripts/classification/output_files/annotationANDlncRNAs.gtf >> T2G.tsv

# Extract ribosomic genes from gtf to filter them from quantification
grep 'gene_biotype "rRNA"' 00_scripts/classification/output_files/annotationANDlncRNAs.gtf | \
	awk -F 'gene_id "' '{ if (NF>1) {split($2, a, "\""); print a[1]}}' > ref/ribo.txt
