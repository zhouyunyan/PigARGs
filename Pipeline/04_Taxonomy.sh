#!/bin/bash

# Kill script if any commands fail
set -e
echo "Job Start at `date`"

mkdir 04_Taxonomy

##############Part4 Taxonomic annotation###########################
diamond=~/bin/diamond
basta=~/basta
Database=~/Database
GeneCatalog=./03_Gene_Catalog
Taxonomy=./04_Taxonomy


#make nr database
diamond makedb --in ${Database}/nr/nr -d ${Database}/nr/nr

#aligning protein sequence of gene catalog to the nr database
${diamond} blastp -q ${GeneCatalog}/02_cdhit_cluster/PIGC_prot.fa -d ${Database}/nr/nr.dmnd -t tmp -p 8 -e 1e-5 -k 50 --id 30 --sensitive -o ${Taxonomy}/PIGC_prot.fa.diamond2nr    

#taxonomic classification based on the LCA algorithms
uniprot=~/Database/basta/taxonomy/uniprot_mapping.db
${basta} sequence -l 25 -i 50 -e 0.00001 -m 3 -b 1 -p 60 ${Taxonomy}/PIGC_prot.fa ${Taxonomy}/PIGC_prot.fa.dmnd.lca.out ${Database}/basta/uniprot

#get time end the job
echo "Job finished at:" `date`
