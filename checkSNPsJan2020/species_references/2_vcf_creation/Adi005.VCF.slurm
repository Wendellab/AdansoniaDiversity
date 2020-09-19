#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70 
#SBATCH --mem=500G 
#SBATCH --time=3-02:30:00
#SBATCH --output=job%j.Adi005VCF.%N.out 
#SBATCH --error=job%j.Adi005VCF.%N.err
#SBATCH --job-name="Adi005VCF" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 

name=Adi005
ref=$name.reference.fasta
thr=70

module load sentieon-genomics/201808.01-opfuvzr

time ( sentieon driver -t $thr -r $ref -i $name.realign.bam  --algo Genotyper $name.VCF  ) &> $name.vcf.timelog

