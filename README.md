# biofilm-data-analysis

The repository contains data analysis scripts used in the paper [Embryo-like features in developing Bacillus subtilis biofilms](https://www.biorxiv.org/content/10.1101/2020.03.09.983718v1).

##Requirements
To run the scripts you need to install Python 3.6 or later and R V3.4.2 or later.
Available at:
- http://www.python.org/getit/
-  https://www.r-project.org/

## Transcriptome data analysis
The [Bs3610_biofilm_transcriptome_analysis.R](scripts/Bs3610_biofilm_transcriptome_analysis.R) script makes raw counts for biofilm transcriptome timeseries expression dataset, transcriptome PCA, pairwise differential expression and likelihood ratio test.

The input files are:
- bam files (obtained from original raw sequences deposited at NCBI GEO (accession number GSE141305)
- Bs3610 gff file from NCBI GenBank (see [Bs3610.gff](example_data/trans_analysis/Bs3610.gff))
- Bs3610 chromosome size file (see [Bs3610.chrsizes](example_data/trans_analysis/Bs3610.chrsizes))
- file containing timepoint and replicate names (see [Sample_Name_biofilm.txt](example_data/trans_analysis/Sample_Name_biofilm.txt))

## Fraction of transcripts calculation
The [raw_by_length.R](scripts/raw_by_length.R) script makes file with raw counts divided by their gene lengths for genes that passed the phylostratigraphic procedure. 

### Example
The input files are:
- raw counts table obtained from DESeq2 [Bs3610_biofilm_transcriptome_analysis.R](scripts/Bs3610_biofilm_transcriptome_analysis.R) (see [Bs3610_biofilm_raw_counts.txt](example_data/raw_by_length/Bs3610_biofilm_raw_counts.txt))
- file containing protein_id, locus_tag and gene_length obtained from NCBI GenBank Bs3610 feature table (see [protein_id_locus_tag_gene_length.txt](example_data/raw_by_length/protein_id_locus_tag_gene_length.txt))
- file containing protein_id and phylostratum obtained from phylostratigraphic procedure ([phyl_tree.py](scripts/phyl_tree.py)) (see [map_full_blast.txt](example_data/phyl/results/maps/map_full_blast.txt))

```console
# help
Rscript raw_by_length.R -h

# run script
Rscript raw_by_length.R ../example_data/raw_by_length/Bs3610_biofilm_raw_counts.txt ../example_data/raw_by_length/protein_id_locus_tag_gene_length.txt ../example_data/phyl/results/map/map_full_blast.txt ../example_data/raw_by_length/results/biofilm_raw_by_length.txt
```

## Phylogenetic mapping
The [phyl_tree.py](scripts/phyl_tree.py) module implements methods for mapping each focal species gene to a phylogenetic tree. By default the mapping is performed with the phylostratigraphy method (blast input), a variation of the method that is based on gene clusters can be used with flag **--cluster**.
When using the default mapping option the module can perform the analysis of gaps in the phylogenetic tree mappings (flag **--hgt**). The gaps can be used as indications of horizontal gene transfer, the genes can be filttered or their phylostrata can be corrected.

Run module with -h for detailed information on functionalities and parameters.

### Example
The input files are:
- tree nodes file: [example_data/phyl/nodes_Bacillus_biofilm.dmp](example_data/phyl/nodes_Bacillus_biofilm.dmp)
- tree names file: [example_data/phyl/names_Bacillus_biofilm.dmp]([example_data/phyl/names_Bacillus_biofilm.dmp)
- blast results file: [example_data/phyl/blast_example](example_data/phyl/blast_example)

```console
# help
python3 phyl_tree.py -h

# run basic mapping
python3 phyl_tree.py ../example_data/phyl/nodes_Bacillus_biofilm.dmp ../example_data/phyl/names_Bacillus_biofilm.dmp 535026 ../example_data/phyl/biofilm_blast_output ../example_data/phyl/results/maps/map_full_small.txt

# run mapping gap counts, export filtered and corrected maps
python3 phyl_tree.py ../example_data/phyl/nodes_Bacillus_biofilm.dmp ../example_data/phyl/names_Bacillus_biofilm.dmp 535026 ../example_data/phyl/biofilm_blast_output ../example_data/phyl/results/ --ps_merge ../example_data/phyl/ps_merged.txt --hgt
```
Basic mapping result file:
- mappings of genes phylostratigraphy method (blast input): [example_data/phyl/results/maps/map_full_blast.txt](example_data/phyl/results/maps/map_full_blast.txt) : results for full blast input (not example)

Horizontal gene transfer analysis (flag **--hgt**) results files ([example_data/phyl/results/hgt/](example_data/phyl/results/hgt/)):
- [ps_species_count.xlsx](example_data/phyl/results/hgt/ps_species_count.xlsx) : number of species per phylostratum (ps)
- [gene_hit_species_count.xlsx]((example_data/phyl/results/hgt/gene_hit_species_count.xlsx) : number of species hits per gene per phylostratum (here gaps can be seen)
- [gene_hit_species_percent.xlsx](example_data/phyl/results/hgt/gene_hit_species_percent.xlsx) : percent of species hits (number_species_hits / total_number_species_ps)
- [gene_hit_species_stats.xlsx](example_data/phyl/results/hgt/gene_hit_species_stats.xlsx) : percent hits per ps, original ps, hgt correction ps, hgt iterative correction ps, total species hit perc (from original ps to focal species), gap size (used in correction, flag **--hgt_gap**)
- [map_filter_2.txt](example_data/phyl/results/hgt/[map_filter_2.txt) : map with removed genes with one or more gaps grater or equal to zero filter size (flag **--hgt_zero_filter**, default 2)
- [map_hgt_correction_2.txt](example_data/phyl/results/hgt/map_hgt_correction_2.txt) : map with genes ps corrected for a single gap (only first) grater or equal to hgt gap size (flag **--hgt_gap**, default 2)
- [map_hgt_correction_iter_2.txt](example_data/phyl/results/hgt/map_hgt_correction_iter_2.txt) : map with genes ps iteratively corrected for gaps grater or equal to hgt gap size (flag **--hgt_gap**, default 2)

## Processing of time series expression data
The module [process_expression_time_series.py](scripts/process_expression_time_series.py) processes gene expression time series values with replicates.
Depending on the given arguments different steps can be included (default 1 and 2).
Process steps:
1. Calculate frequencies per replicate.
2. Resolve replicates to one value, replicate median is calculated from non zero replicate values.
3. Zero filtering, genes with zero values for two or more points in time series are removed.
4. Interpolation of zero values with neighbours mean or with neighbour value if zero is located at the start or end of the series (connected zero values are not allowed).
5. Median log<sub>2</sub> transformation per gene time series. All values in gene expression profiles are divided by the gene expression profile median value, after which the log<sub>2</sub> transformation is performed.

Run module with -h for detailed information on parameters.

### Example
The input file is a tab delimited text file with expression values, time point names must be separated with "_" from replicate numbers ([example_data/time_series_expression/trans_replicates_norm_length.txt](example_data/time_series_expression/trans_replicates_norm_length.txt))
```console
# help
python3 process_expression_time_series.py -h

# run with default params (used for TAI, PAI calculation)
python3 process_expression_time_series.py  ../example_data/time_series_expression/trans_replicates_norm_length.txt ../example_data/time_series_expression/trans_TAI_res.txt 

# run with all process steps (used for profile plotting)
python3 process_expression_time_series.py ../example_data/time_series_expression/trans_replicates_norm_length.txt ../example_data/time_series_expression/trans_profiles_res.txt --transform --double_zero_remove --interpolate
```

## Functional analysis
The module [functional_analysis.py](scripts/functional_analysis.py) module performs enrichment analysis of gene functions per phase.
The enrichment analysis is performed with a hypergeometric test, the multiple testing correction is done with the Benjamini Hochberg method.
The enrichment test is based on annotations which are combinations of genes and corresponding functions.
The hypergeometric test for a function *X*  in a phase *P*  is performed on counts:
- number of *X* function annotations in phase *P* (quant)
- number of annotations in phase *P* (sample)
- number of *X* function annotations in all phases (hit)
- total number of annotations (total)

The module implements the analysis for different ontology depths and expression cutoff values (see config file).

### Example
The inputs are gene category, expression and gene feature files (see [example_data/functional_analysis](example_data/functional_analysis)).
The arguments (file paths, cutoffs, depths...) are given through a config file, see [example_data/functional_analysis/config.txt](example_data/functional_analysis/config.txt) for an example and instructions.

```console
# help
python3 functional_analysis.py -h

# run
python3 functional_analysis.py ../example_data/functional_analysis/config.txt
```

#### Result files
The result files are organised in folders by expression cutoff values and ontology depths. For example, the results for cutoff value 1 and dept 3 will be located in folder : [example_data/functional_analysis/results/cutoff_1_0/depth_3](example_data/functional_analysis/results/cutoff_1_0/depth_3/) : 
- [hyper_cutoff_1_0_depth_3.txt](example_data/functional_analysis/results/cutoff_1_0/depth_3/hyper_cutoff_1_0_depth_3.txt)
Main results file, contains quant, sample, hit, total, p_value, p_adj, nad log_odds for every function in every phase.
- [hyper_filtered_cutoff_1_0_depth_3.txt](example_data/functional_analysis/results/cutoff_1_0/depth_3/hyper_filtered_cutoff_1_0_depth_3.txt)
Contains only significant (sig_p param in config file) results from main results file.
- [p_adj_values_cutoff_1_0_depth_3.xlsx](example_data/functional_analysis/results/cutoff_1_0/depth_3/p_adj_values_cutoff_1_0_depth_3.xlsx)
Adjusted p-values of a hypergeometric test in a table format (used for heatmap plotting).
- [log_odds_cutoff_1_0_depth_3.xlsx](example_data/functional_analysis/results/cutoff_1_0/depth_3/log_odds_cutoff_1_0_depth_3.xlsx)
Log odds values in a table format.
- [annotations_all.txt](example_data/functional_analysis/results/cutoff_1_0/depth_3/annotations_all.txt)
All annotations, an annotation is a combination of gene-id and its function.
- [annotations_phase_cutoff_1_0_depth_3.txt](example_data/functional_analysis/results/cutoff_1_0/depth_3/annotations_phase_cutoff_1_0_depth_3.txt)
Annotations per phase.
- [genes_functions_cutoff_1_0_depth_3.xlsx](example_data/functional_analysis/results/cutoff_1_0/depth_3/genes_functions_cutoff_1_0_depth_3.xlsx)
The file contains gene-ids and gene names for all significant functions per phase.

## Transcriptome evolutionary indices
The [trans_evo_ind.R](scripts/trans_evo_ind.R) script calculates and plots transcriptome evolutionary indices Transcriptome Age Index (TAI), Transcriptome nonsynonymous divergence index (TdNI), Transcriptome synonymous divergence index (TdSI) and Transcriptome Codon Bias Index (TCBI). 

The input files are:
- file with locus_tag and fraction of transcripts and replicates resolved (obtained by [raw_by_length.R](scripts/raw_by_length.R) and [process_expression_time_series.py](scripts/process_expression_time_series.py) steps 1 and 2) for calculating all transcriptome evolutionary indices (see [biofilm_raw_by_length.txt](example_data/trans_evo_ind/biofilm_raw_by_length.txt))
- file containing protein_id, locus_tag and gene_length obtained from NCBI GenBank Bs3610 feature table (see [protein_id_locus_tag_gene_length.txt](example_data/raw_by_length/protein_id_locus_tag_gene_length.txt)) for matching protein_id with locus_tag
- file containing protein_id and phylostratum obtained from phylostratigraphic procedure ([phyl_tree.py](scripts/phyl_tree.py)) (see [map_full_blast.txt](scripts/example_data/phyl/results/maps/map_full_blast.txt)) for calculating TAI
- fasta file of Bs3610 and BlDSM13 protein sequences obtained from NCBI GenBank (see [Bs3610_protein.faa](example_data/trans_evo_ind/Bs3610_protein.faa) and [BlDSM13_protein.faa](example_data/trans_evo_ind/BlDSM13_protein.faa)) for calculating dNdS, TdNI and TdSI
- fasta files of Bs3610 cds and ribosomal sequences obtained from NCBI GenBank (see [Bs3610_cds_from_genomic.fna](example_data/trans_evo_ind/Bs3610_cds_from_genomic.fna) and [Bs3610_ribosomal.fasta](example_data/trans_evo_ind/Bs3610_ribosomal.fasta)) for calculating codon usage bias and TCBI

The output files you will get are:
- TAI plot
- TdNI plot
- TdSI plot
- TCBI plot


## Team

* **Momir Futo**
* **Luka Opašić**
* **Sara Koska**
* **Nina Čorak**
* **Tin Široki**
* **Vaishnavi Ravikumar**
* **Annika Thorsell**
* **Domagoj Kifer**
* **Mirjana Domazet-Lošo**
* **Kristian Vlahoviček**
* **Ivan Mijaković**
* **Tomislav Domazet-Lošo**

## Acknowledgments
This work was supported by City of Zagreb Grant, Croatian Science Foundation under the project IP-2016-06-5924, Adris Foundation Grant and European Regional Development Fund Grants KK01.1.1.01.0008 (CERRM) and KK.01.1.1.01.0009 (DATACROSS).