# Pipeline for metagenomic analysis

This directory contains scripts related to the manuscript "Extensive metagenomic analysis of the porcine gut resistome to identify indicators reflecting
antimicrobial resistance". 

Before running, you must ensure that all required softwares and databases are installed successfully. 

## INSTALLATION

Create two directories "bin" and "Database" in user home directory. 

### Software installation

The installation method refer to the manual of each software. The name, version and availablity of the software are as follows:  

|Software|Availability|
|:-----|:---------|
|fastp (v0.19.4)|https://github.com/OpenGene/fastp|
|bwa (v0.7.17-r1188)|https://github.com/lh3/bwa|
|Samtools (v1.10)|https://github.com/samtools/samtools/releases/|
|bedtools (v2.28.0)|https://bedtools.readthedocs.io/en/latest/|
|MEGAHIT (v1.1.3)|https://github.com/voutcn/megahit|
|Bowtie 2 (v2.3.4.1)|https://anaconda.org/bioconda/bowtie2|
|MetaGeneMark (v3.38)|http://exon.gatech.edu/heuristic_gmhmmp.cgi|
|CD-HIT (v4.7)|https://github.com/weizhongli/cdhit/releases|
|featurecount (v2.0.1)|http://bioinf.wehi.edu.au/featureCounts/|
|diamond (v0.9.24)|http://www.diamondsearch.org/index.php|
|BASTA (v1.3.2.3)|https://github.com/timkahlke/BASTA|
|RGI (v5.1.1)|https://card.mcmaster.ca|

Note: Make all needed command of software availabled in the "~/bin" directory or in system environment variables.The version is only the version used in the paper and does not have to be the same. 

### Database installation

All databases are stored in the "~/Datebase" directory. 

The name,description and availavlity of the database are as follows: 

|Database|Version/release date|Description|Availability|
|:-------|:-------------------|:----------|:-----------|
|Pig (Sscrofa11.1)|Sscrofa11.1|Pig reference genome|http://asia.ensembl.org/Sus_scrofa/Info/Index|
|NCBI NR|version 2019_04|protein database|ftp://ftp.ncbi.nlm.nih.gov/blast/db|
|CARD|v3.1.0|Antibiotic Resistance genes annotation|https://github.com/arpcard/rgi#install-dependencies|

Note: The version are only the version used in the paper,most of database are constantly updated.

## OVERVIEW OF PIPELINE

The scripts of metagenomic analysis are placed in "[Pipeline](https://github.com/zhouyunyan/PigARGs/tree/main/Pipeline)" directory. The scripts of statistical analysis and visualization are placed in "[Scripts](https://github.com/zhouyunyan/PigARGs/tree/main/Scripts)" directory.

### Metagenomic analysis

#### Part1: 01_data_preprocessing.sh

Metagemonic data pre-processing: read trimming and host (pig) read removal,generating high-quality sequence. 

#### Part2: 02_Assembly.sh

Metagenomic assembly: Assemble short reads into long contigs.

A total of four scripts in this modules, including gene prediction, taxonomy annotation, function annotation and abundance estimation.

#### Part3: 03_Gene_Catalog.sh 

This part contains steps of gene prediction, filtration of incomplete genes, integration of gene catalog and gene dereplications.

#### Part4: 04_Taxonomy.sh 

The protein sequence of genes were aligned to NCBI NR database, and the taxonomic classification were determined based on the last (lowest) common ancestor algorithms.

#### Part5: 05_Function.sh

The antibiotic resistance genes were identified by alignment against the Comprehensive Antibiotic Resistance Database (CARD).

#### Part6: 06_Abundance.sh

Gene abundance were caculated by aligning clean reads of each sample to the gene catalog to obtained the counts of mapped reads, and  normalized to read count fragments per kilobase million (FPKM). The abundance of function items were performed were calculated by adding the abundances of all its members falling within each category with R scripts. 

### Statistical analysis and visualization

Statistical analysis and visualization were handled by scripting with R program. These scripts were placed in "[Scripts](https://github.com/zhouyunyan/PigARGs/tree/main/Scripts)" directory. All related input data for statistical analysis and visualization are in "[Pre-processed_Files](https://github.com/zhouyunyan/PigARGs/tree/main/Pre-processed_Files)" directory.
