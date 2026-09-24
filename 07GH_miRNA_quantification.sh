#!/bin/bash

# Choose base directory
base="/media/labo/Expansion/mir_lymphnode"
cd $base

# FASTQC: We perform the quality control of the raw data to asses the general status of the reads.
fastqc rawfq/*.fastq.gz
multiqc rawfq/*.zip
mv multiqc_data multiqc_report.html log/

# TRIMMOMATIC: Here, we remove the short reads (less than 36 bp) and the bases or reads with low quality with.
cd rawfq
for f1 in *fastq.gz 
do
	f2=${f1%%fastq.gz}"trim_fastq.gz"
	f3=${f1%%fastq.gz}"_summ.log"
cd $base
trimmomatic SE  -phred33 rawfq/$f1 trimfq/$f2 \
	ILLUMINACLIP:ref/adapters.fa2:30:7 LEADING:20 TRAILING:20 \
	SLIDINGWINDOW:5:20 MINLEN:16 &> trimfq/log/trim_log/$f3
echo "${f1%%.fastq.gz} done, see log file: $f3"
done

# FASTQC: We check that low quality bases/reads were removed correctly.
echo "Performing QC of trimmed reads"
fastqc trimfq/*.fastq.gz
multiqc trimfq/*.zip
mv multiqc_data multiqc_report.html log/

# Download the sheep reference genome in fasta format and convert it to bowtie format
gzip -d 00_ref/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa.gz
bowtie-build -f 00_ref/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa oar_index
mkdir 00_references/oar_index
mv *.ebwt 00_references/oar_index/

# Map reads against sheep genome index
for f1 in *trimfq
do
	f2=${f1%%trim_fastq}"collapsed.fa"
	f3=${f1%%trim_fastq}"vs_oar.arf"
   echo "mapping $f1 sample and generating $f2 and $f3"
   mapper.pl trimfq/$f1 -e -h -j -m -p $base/00_references/oar_index -q -s mapping/$f2 -t mapping/$f3
done

# Set a sample identifier to each read in the fasta file and merge all fasta files
cd mapping
while IFS="," read -r col1 col2 col3
do
	fasta=$col1".fa"
	echo $fasta
	sed -i 's/seq/'$col2'/g' $fasta
done < sample_info.csv
cat *collapsed.fa > merged.fa
cd $base

# Quantify miRNAs
quantifier.pl -p 00_references/premirna_seq.fa \
	-m 00_references/mature_seq.fa \
	-r mapping/merged.fa \
	-s 00_references/star_seq.fa \
	-y now -j > mapp_stats.out
