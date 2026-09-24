# sheep_lymph_nodes_transcriptome

This repository contains the scripts used in the manuscript: "Integrated mRNA and ncRNA transcriptomic analysis of lymph nodes in sheep after repetitive vaccination"

## Abbreviation:
* RNA-seq: RNA sequencing
* lncRNA: Long non-coding RNA
* miRNA: microRNA
* DE: Differentially expressed

## Scripts:
* [01GH_Preprocessing_and_lncRNA_prediction.sh](/01GH_Preprocessing_and_lncRNA_prediction.sh): This script performs quality control and preprocessing of the RNA-seq data. Also is used for genome alignement, potential lncRNA selection and coding potential assesment.
  
* [02GH_lncRNA_intersection.R](/02GH_lncRNA_intersection.R): This script performs the intersection of the three coding potential assesment tools and exports a final lncrna candidate list and a Venn diagram.
  
* [03GH_quantification.sh](/03GH_quantification.sh): This script performs the quantification of transcripts using the annotated transcriptome and the discovered novel lncRNAs
  
* [04GH_Expression_analyses.R](/04GH_Expression_analyses.R): This script filter the lowly expressed genes, correct the counts and performs the differential expression analysis.
  
* [05GH_Functional_annotation.R](/05GH_Functional_annotation.R): This script performs the functional annotation of the differentially expressed genes and performs a TreeMap analysis.
  
* [06GH_coexpression.R](/06GH_coexpression.R): This script detect co-expressed gene modules and their hubgenes. Also export node and edge tables to construct a co-expression network in Cytoscape.
  
* [07GH_miRNA_quantification.sh](/07GH_miRNA_quantification.sh): This script performs quality control, preprocessing, alignement and quantification of the miRNA-seq data.
  
* [08GH_mirna_expression_analysis.R](/08GH_mirna_expression_analysis.R): This script filter the lowly expressed miRNAs and performs the differential expression analysis.
  
* [09GH_target_prediction.R](/09GH_target_prediction.R): This script performs the intersection of the miRNA target prediction tools and performs the functional annotation of the targets of each DE miRNAs
  
* [10GH_mirna_target_corr.R](/10GH_mirna_target_corr.R):
