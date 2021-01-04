#!/bin/bash

# Kill script if any commands fail
set -e
echo "Job Start at `date`"

mkdir 05_Function

##############Part5 Functional annotation###########################
cdhit=./03_Gene_Catalog/02_cdhit_cluster/
Database=~/Database
faa=PIGC_prot.fa
Function=./05_Function

rgi=~/bin/rgi

cd ${Function}
mkdir CARD

#CARD
$rgi main -i ${cdhit}/$faa -o CARD/$faa.card -n 6 --debug -t protein -a DIAMOND --clean

#get time end the job
echo "Job finished at:" `date`