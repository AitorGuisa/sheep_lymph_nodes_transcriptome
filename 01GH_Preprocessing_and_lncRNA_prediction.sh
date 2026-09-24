#!/bin/bash

# Choose base directory
base="/media/labo/Expansion/lymph_node_total"
cd $base

# FASTQC: We perform the quality control of the raw data to asses the general status of the reads
fastqc rawfq/*.fastq.gz
multiqc rawfq/*.zip
mv multiqc_data multiqc_report.html log/

# TRIMMOMATIC: Here, we remove the short reads (less than 36 bp) and the bases or reads with low quality
cd rawfq
for f1 in *R1_001.fastq.gz
do
     f2=${f1%%R1_001.fastq.gz}"R2_001.fastq.gz"
     f3=${f1%%R1_001.fastq.gz}"_1_p.fastq.gz"
     f4=${f1%%R1_001.fastq.gz}"_1_up.fastq.gz"
     f5=${f1%%R1_001.fastq.gz}"_2_p.fastq.gz"
     f6=${f1%%R1_001.fastq.gz}"_2_up.fastq.gz"
     f7=${f1%%R1_001.fastq.gz}"_summ.log"
cd $base
trimmomatic PE -threads 8 -phred33 \
	rawfq/$f1 rawfq/$f2 trimfq/$f3 trimfq/$f4 trimfq/$f5 trimfq/$f6 \
	ILLUMINACLIP:ref/adapters.fa:2:30:10:1:"true" LEADING:20 TRAILING:20 \
	SLIDINGWINDOW:5:20 MINLEN:36 &> trimfq/log/trim_log/$f7
echo "${f1%%_1.fastq.gz} done, see log file: $f7"
done

# FASTQC: We check that low quality bases/reads were removed correctly
echo "Performing QC of trimmed reads"
fastqc trimfq/*.fastq.gz
multiqc trimfq/*.zip
mv multiqc_data multiqc_report.html log/

# BBTOOLS: Here we remove the ribosomic sequences with bbduk
while read f1
do   
     f2=$f1"_1_p.fastq.gz"
     f3=$f1"_2_p.fastq.gz"
     f4=$f1"_stats.txt"
bbduk.sh -Xmx60000m in=trimfq/$f2 in2=trimfq/$f3 out=rless/$f2 out2=rless/$f3 ref=ref/ribokmers.fa k=31 hdist=1 stats=log/bbduk_log/$f4
done < script/names.txt

# FASTQC: We check the general status of the samples after removing ribosomic sequences
echo "Performing QC of trimmed reads"
fastqc rless/*.fastq.gz
multiqc rless/*.zip
mv multiqc_data multiqc_report.html log/

# STAR: We generate the indexes for the alignment and we introduce the index and the trimmed samples in STAR to align the transcripts
STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir ref/sindex/ \
--genomeFastaFiles ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna \
--sjdbGTFfile ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.gtf \
--sjdbOverhang 107 &> log/sindex.log
mkdir /tmp/temp_star
chmod 777 /tmp/temp_star

while read f1
	do
	f2=$f1"_1_p.fastq.gz"
	f3=$f1"_2_p.fastq.gz"
	f4=$f1"."
	echo "Mapping: $f1"		
	STAR --runThreadN 8 \
	--genomeDir ref/sindex/ \
	--readFilesIn rless/$f2 rless/$f3 \
	--readFilesCommand zcat \
	--outFileNamePrefix mapping/$f4 \
	-outSAMtype BAM SortedByCoordinate \
	--outMultimapperOrder Random \
	--outSAMattrIHstart 0 \
	--twopassMode Basic \
	--outTmpDir /tmp/temp_star/temp &> log/mapping/$f1".log"
	echo "$f1 done, see log file $f1 .Log.final.out"
done < names.txt

# MULTIQC and QUALIMAP: we perform quality controls of the genome alignement
multiqc mapping/*Aligned.sortedByCoord.out.bam
mv multiqc_data multiqc_report.html log/mapping/
cd mapping
for f1 in *Aligned.sortedByCoord.out.bam
do
	f2=${f1%%_TOTAL_MERGED_S1_L001.Aligned.sortedByCoord.out.bam}
	cd $base
	mkdir -p log/mapping/mapqc/$f2
	qualimap rnaseq -bam mapping/$f1 -gtf ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.gtf -outformat HTML -outdir log/mapping/mapqc/$f2 -p strand-specific-reverse -pe --java-mem-size=12G    
	echo "$f2 done"
done

# STRINGTIE: we assemble de sheep transcriptome de novo to detect new potential lncRNA transcripts
cd $base/mapping
# we must remove "gene" entries as stringtie and gffcompare don't accept empty transcript id values
awk '$3 !="gene"' ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.gtf > ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic_nogenes.gtf 

for f1 in *Aligned.sortedByCoord.out.bam
do
	f2=${f1%%_TOTAL_MERGED_S1_L001.Aligned.sortedByCoord.out.bam}".gtf"
	cd $base
	echo "Assembly of sample $f2 started."
	stringtie mapping/$f1 -o assembly/$f2 -p 8 -G ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic_nogenes.gtf --rf
	echo "Assembly of sample $f2 ended"
done 

# merge all assemblies into a non-redundant annotation
cd $base/assembly
ls *.gtf  > gtf_list.txt
stringtie --merge -G /media/labo/Expansion/lymph_node_total/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic_nogenes.gtf -o merged.gtf gtf_list.txt

# GFFCOMPARE: Compare the merged annotation with the available annotation 
gffcompare -r "/media/labo/Expansion/lymph_node_total/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic_nogenes.gtf" -V "/media/labo/Expansion/lymph_node_total/assembly/merged.gtf" -o comparison

# lncRNA Characterization: Select candidate lncRNAs as those classified as "u" (unknown intergenic), "i" (intronic),"x" (antisense) and "o" (other same strand overlap)
awk -F "\t" 'BEGIN{OFS="\t"} $3=="u" || $3=="i" || $3=="x" || $3=="o" {print $0}' comparison.merged.gtf.tmap > lncRNA_candidate.tmap

# Filter out those transcripts whose length is <200nt if they are multiexonic or <2000nt if they are single-exon
awk -F "\t" 'BEGIN{OFS="\t"} {if( $6 >= 2 && $10 >= 200 ) {print $0} else if( $6 == 1 && $10 >= 2000 ) {print $0}}' lncRNA_candidate.tmap > lncRNA_candidate2.tmap && mv -f lncRNA_candidate2.tmap lncRNA_candidate.tmap

# Check their coding potential and if they contain protein domains with multiple tools, and select those that pass all the filtering criteria. First we have to extract the nt sequence of the candidate lncRNAs
awk -F "\t" '{print$5}' lncRNA_candidate.tmap > lncRNA_list.txt
awk -F "\t" 'FNR==NR {myarray[$1]++; next} {aux=$9; gsub(/^transcript_id "/,"",aux); gsub(/".+/,"",aux); if( myarray[aux] ) {print $0}}' lncRNA_list.txt comparison.annotated.gtf > lncRNA_candidate.gtf
gffread -w lncRNA_candidate.fa -g "/media/labo/Expansion/lymph_node_total/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna" lncRNA_candidate.gtf
mv lncRNA_candidate.tmap lncRNA_candidate.gtf lncRNA_list.txt lncRNA_candidate.fa ../lnc

# CPC2: Check coding potential 
/home/labo/soft/CPC2_standalone-1.0.1/bin/CPC2.py -i lncRNA_candidate.fa -o CPC2_results
sed '1d' CPC2_results.txt | awk -F "\t" '$8=="noncoding" {print $1"\tncRNA"} $8=="coding" {print $1"\tprotein_coding"}' > CPC2_results_summary.txt

# CPAT: Check coding potential.Apply model to detected lncRNAs, using a cow model
modelr="/media/labo/Datuak/Encefalo/cpat/cowmodel/"
modelhex="/media/labo/Datuak/Encefalo/cpat/cowmodel/"
mkdir cpat
mv lncRNA_candidate.fa cpat/
cd cpat
cpat.py -g lncRNA_candidate.fa -d $modelr/cow.logit.RData -x $modelhex/cow_hexamer.tsv -o CPAT_results
awk -F "\t" '{print $1"\tncRNA"}' CPAT_results.no_ORF.txt > CPAT_results.txt

#CPAT is very strict with genes that do not have an ORF, so we select the genes with a probability below 0.349
sed '1d' CPAT_results.ORF_prob.best.tsv | awk -F "\t" '$11>=0.349 {print $1"\tprotein_coding"} $11<0.349 {print $1"\tncRNA"}' >> CPAT_results.txt
mv CPAT_results.txt ..
cd ..

# HMMER: Check for protein domains
# First translate all lncRNA candidates into their 3 possible Open Reading Frames (ORFs) with transeq
transeq -sequence lncRNA_candidate.fa -frame F -outseq transeq_results.fa
# Then, download the pfam-A database and check againts it if there is any protein domain with hmmscan (from HMMER)
hmmpress Pfam-A.hmm
hmmscan --tblout lncRNA_hmmscan.persequence.txt  --noali "/media/labo/Expansion/lymph_node_total/ref/pfam/Pfam-A.hmm" transeq_results.fa
# Select those with a domain hit with E-value<1e-5 in full sequence ($5) and best 1 domain ($8)
sed '1d' lncRNA_hmmscan.persequence.txt | awk '$5<1e-5 && $8<1e-5 {gsub(/_[0-9]/,"",$3); print $3}' | sort -u > list.txt
awk -F "\t" '{print $1}' CPAT_results.txt > tmp.txt
awk -F "\t" 'FNR==NR {myarray[$1]++; next} myarray[$1] {print $1"\tprotein_coding"; next} {print $1"\tncRNA"}' list.txt tmp.txt > hmmscan_results.txt
rm list.txt tmp.txt
