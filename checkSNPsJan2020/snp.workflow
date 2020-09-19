#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 168:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 100   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "AdanAlleles"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-Adansonia.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

ml samtools
ml parallel 

mkdir checkSNPsJan2020/
cd checkSNPsJan2020/

for a in ../bam/*.bam; do basename -s .sort.resort.bam $a >> bam.list; done

parallel -j 10 'samtools bam2fq ../bam/{}.sort.resort.bam > {}.fastq' :::: bam.list 
mkdir fasta && cp ../output/fasta/*.nonredundant.fasta fasta
cat fasta/* | grep ">" | sed '/Gorai/d' | sort | uniq > alleles.list

cd fasta
for a in *Gorai*.fasta; do sed "s/>/>${a%.nonredundant.fasta}_/g" $a >> all.sequences.fasta; done

mkdir ../species_references/
cat ../bam.list | while read line; do ./../seqkit grep -p .*$line -r all.sequences.fasta > ../species_references/$line.reference.fasta; done

cd ../species_references/
ml bwa

ls *.fasta | parallel -j 20 bwa index {}

# after this run (in species_references):
# slurm files in 1_sentieon  
# slurm files in 2_vcf_creation
# extract.info.bash
# make.table.slurm