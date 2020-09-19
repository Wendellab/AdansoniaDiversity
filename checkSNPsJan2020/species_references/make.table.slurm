#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70 
#SBATCH --mem=500G 
#SBATCH --time=3-02:30:00
#SBATCH --output=job.table.%J.out 
#SBATCH --error=job.table.%J.err
#SBATCH --job-name="table" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 

module load parallel
module load r/3.6.0-py2-fupx2uq

# uses two files generated from VCF to make a table of SNP positions, genotypes, number of reads, gene length, and %gene with variants
# snplist file : for a in *.VCF; do cut -f1,2,10 $a | sed '/#/d' | cut -f1,2 -d ':' | sed 's/[:,]/\t/g' > ${a%.VCF}.snplist; done
# allele.length file : for a in *.fai; do cut -f1,2 $a > ${a%%.*}.allele.length; done

ls *.VCF | parallel -j 10 'Rscript make.variant.table.R {.}.snplist {.}.allele.length {.}.tbl'


