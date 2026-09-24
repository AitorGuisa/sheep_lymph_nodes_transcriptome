# sheep_lymph_nodes_transcriptome

This repository contains the scripts used in the manuscript: 

## Scripts
* [01GH_Preprocessing_and_lncRNA_prediction.sh](/01GH_Preprocessing_and_lncRNA_prediction.sh):
  
* [02GH_lncRNA_intersection.R](/02GH_lncRNA_intersection.R):
  
* [03GH_quantification.sh](/03GH_quantification.sh):
  
* [04GH_Expression_analyses.R](/04GH_Expression_analyses.R):
  
* [05GH_Functional_annotation.R](/05GH_Functional_annotation.R):
  
* [06GH_coexpression.R](/06GH_coexpression.R):
  
* [07GH_miRNA_quantification.sh](/07GH_miRNA_quantification.sh):
  
* [08GH_mirna_expression_analysis.R](/08GH_mirna_expression_analysis.R):
  
* [09GH_target_prediction.R](/09GH_target_prediction.R):
  
* [10GH_mirna_target_corr.R](/10GH_mirna_target_corr.R):
  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXAMPLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* [mirdeep2-core-command.sh](/mirdeep2-core-command.sh): Code used to run the preprocessing of the samples, genome mapping and the core miRDeep2 algorithm.

* [mirdeep2-quantifier.sh](/mirdeep2-quantifier.sh): Code used to run the miRDeep2 quantifier algorithm.

* [Statistic_analysis.R](/Statistic_analysis.R): Code used for the analysis of miRNA expression and tissue specificity.

* [novel_mirnas.R](/novel_mirnas.R): Code used to give correct "3p" and "5p" names to the filtered miRNAs and to prepare the mature and pre-miRNA fasta files for quantification.  

* [mirna_blast.py](/mirna_blast.py): Code used for sequence conservation analysis of novel miRNAs, selection of unique miRNA sequences for quantification and search of clusters in genome.

* [expression_plots.py](/expression_plots.py): Code used to plot miRNA expression by conservation status, tissue and specificity.
