# get the info we need for the tables
mkdir dedupBam joblogs logs metrics realignBam references sam slurmFiles sortBam stats
mv *.realign.bam* realignBam
mv *.sort.bam* sortBam
mv *.dedup.bam* dedupBam
mv job.* joblogs
mv *.slurm slurmFiles
mv *.sention.slurm slurmFiles/1_sentieon
mv *.VCF.slurm slurmFiles/2_vcf_creation

mv *.metric* metrics
mv *.score* metrics
mv *.summary metrics
mv *.sam sam

for a in *.VCF; do cut -f1,2,10 $a | sed '/#/d' | cut -f1,2 -d ':' | sed 's/[:,]/\t/g' > ${a%.VCF}.snplist; done

for a in *.fai; do cut -f1,2 $a > ${a%%.*}.allele.length; done

mv *.fasta* references