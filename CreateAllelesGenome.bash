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

# load modules

module load parallel/20170322-36gxsog bwa/0.7.17-zhcbtza
#module load gmap-gsnap/2019-05-12-zjqshxf
# module load bambam/1.3-3u4zcom
# module load bedtools2/2.27.1-s2mtpsu
# module load py-biopython/1.70-py2-hcdisd2
module load raxml/8.2.11-openmpi3-ffoixg2
# module load r/3.5.0-py2-x335hrh
# ml bioawk/1.0-3f2ahpq

# from previously hand-curated data, get all the sequences
cd genes_unique_from_344set/
for a in *.fasta; do sed -i -e "/>/s/$/-${a%.fasta}/" $a; done && cd ..
cat 344.baits.fasta genes_unique_from_344set/*fasta > 344.expanded.fasta
 
grep ">" 344.expanded.fasta | cut -f2 -d "-" | sort | uniq > gene.list
ls *.S.fq.gz | cut -f1 -d '.' > sample.list
 
##### Step 1 : Find reasonably covered targets in genome 2017-04-25 release (http://gigadb.org/dataset/100445)(https://doi.org/10.1093/gigascience/giy051)
 
gmap_build -D gmapdb/ -d Bombax Mumian.Scaffold.fasta
gmap -D gmapdb/ -d Bombax -t 50 -n 0 -f 4 344.expanded.fasta > Bombax.gmap.out  # 359 matches
 
mkdir matches_Bombax
parallel "grep {1} Bombax.gmap.out > matches_Bombax/{1}.match" :::: gene.list 
parallel "cut -f1,4,5 matches_Bombax/{1}.match | sort | uniq > matches_Bombax/{1}.reduced" :::: gene.list
# 359 matches parsed
# 
for a in matches_Bombax/*.reduced; do echo $a `cut -f1 $a | sort | uniq | wc -l` >> bombax.matches; done
sed '/ [2340]/d' bombax.matches > bombax.single.matches

cut -f1 -d ' ' bombax.single.matches | while read line; do awk ' { print $0 "\t" FILENAME } ' $line >> bombax.bed; done
sed -i 's/matches_Bombax\///g' bombax.bed


# #### prepare file to extract gene regions
Rscript prepare.bed.R # makes the bed file and filters for genes <10001 nt

cut -f2,3,4,5 bombax.genes.bed | sed 's/^.*S/S/g' | sed '1d' | sed 's/[.]reduced//g' > bombax.extractGenes.bed

bedtools getfasta -fi Mumian.Scaffold.fasta -bed bombax.extractGenes.bed -split -name > Bombax.genes.fasta

#### map trimmed reads to extracted bombax genes
ls *.S.fq.gz | sed 's/[.]S[.]fq[.]gz//g' > reads.list

gmap_build -D gmapdb/ -d BgenesAll Bombax.genes.fasta

parallel -j 3 "gsnap --gunzip -n 1 -Q -t 50 -D gmapdb/ -d BgenesAll -A sam {}.R1.fq.gz {}.R2.fq.gz > {}.P.sam 2>> {}.P.gsnap.log" :::: reads.list

parallel -j 3 "gsnap --gunzip -n 1 -Q -t 50 -D gmapdb/ -d BgenesAll -A sam {}.S.fq.gz > {}.S.sam 2>> {}.S.gsnap.log" :::: reads.list

ls *.sam | parallel -j 10 'samtools view -@ 2 -u -F 4 {} | samtools sort -n -m 3G -o {.}.sort.bam &> {.}.bam.err' 
for a in *.S.sort.bam; do samtools merge ${a%%.*}.sort.bam $a ${a%%.*}.trim.P.sort.bam; done

rm *.[SP].sort.bam
ls *.sort.bam | parallel -j 15 'samtools sort -@ 2 -o {.}.resort.bam {}'
ls *.resort.bam | parallel -j 30 'samtools index {}'

# use bam2consensus from bambam to get the mapped reads back as a fasta file
# outputs one file per gene
bam2consensus -m 5 -p 4 *.resort.bam

# take bam2consensus output and concatenate that with the original Bombax gene (from the reference)
ls *Gorai*.fasta | while read line; do echo ${line%.fasta} > ids.txt; xargs samtools faidx Bombax.genes.fasta < ids.txt >> $line; done
 
# filter the alignments for sequences with <90% and positions with >30% missing data
ls *Gorai*.fasta | parallel -j 20 'sh filter.alignments {}'
rm *.bak *.tmp.fasta

# remove identical sequences from same accession using bioawk
for a in *.Nfil.fasta; do bioawk -c fastx '{print ">"$name" "$comment"\t"$seq}' $a | sed 's/_[1234]//g' | sort -u | sed 's/\t/\n/g' | awk ' /^>/ && _[$1]++ {$1 = $1 "_" (_[$1] + 1 )}; 1' > ${a%.Nfil.fasta}.nonredundant.fasta; done


# mkdir bam intermediate_fasta log reads
# mv *.fq.gz reads/
# mv *.bam* bam/
# mv *.log log/
# mv *.Nfil*fasta intermediate_fasta/
# mv *.rename*fasta intermediate_fasta/
# mv *[0-9].fasta intermediate_fasta/
# 
# # make basic trees
ls *.nonredundant.fasta | parallel -j 100 'raxmlHPC -m GTRGAMMA -p 25632 -s {} -n {.}.GTRGAMMA.raxml -N autoMRE -x 451842 -f a -k'




