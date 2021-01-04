#!/bin/bash

# Kill script if any commands fail
set -e
echo "Job Start at `date`"

mkdir 03_Gene_Catalog

##############Part3 Construction of the gene catalog###########################
CleanData=./01_CleanData
Assembly=./02_Assembly
GeneCatalog=./03_Gene_Catalog
Scripts=~/Scripts
PGC=./03_Gene_Catalog/PGC
gmhmmp=~/bin/gmhmmp 
MetaGeneMark=~/bin/MetaGeneMark_v1.mod
cdhit=~/bin/cdhit

cd ${GeneCatalog}
mkdir 01_before_cdhit 02_cdhit_cluster 03_Gene_abundance

###gene prediction
cd 01_before_cdhit
#gene prediction
$gmhmmp -a -d -f G -p 1 -m $MetaGeneMark ../../${Assembly}/allSample500.final_contigs.fasta -A allSample500_protein.faa -D allSample500_nucl.fna

cd ../../
faa=${GeneCatalog}/01_before_cdhit/allSample500_protein.faa
cd ${GeneCatalog}/02_cdhit_cluster/
###dereplication at protein level (identity:95%,coverage:90%)
$cdhit -i $faa -o PIGC_prot.fa -c 0.95 -n 5 -G 0 -aS 0.9 -M 20000 -T 8 -g 1 -d 0        

cat PIGC_prot.fa |grep "^>"|awk -F ' ' '{print $1}'|awk -F '>' '{print $2}' >PIGC_geneID.list

#extract the nucleotide sequence corresponding to a protein sequence of PIGC by sequence ID
perl ${Scripts}/extract_fabyid.pl PIGC_cds.fna PIGC_geneID.list ../01_before_cdhit/sample500_BGI287_before_redundancy_nofilter_cds.fna

#return to initial directory
cd ../../

#get time end the job
echo "Job finished at:" `date`

