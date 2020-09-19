#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70 
#SBATCH --mem=500G 
#SBATCH --time=3-02:30:00
#SBATCH --output=job%j.Adi007.%N.out 
#SBATCH --error=job%j.Adi007.%N.err
#SBATCH --job-name="Adi007" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 

module load bwa/0.7.17-zhcbtza
module load samtools

DIR=/work/LAS/jfw-lab/corrinne/Adansonia/AdansoniaAllelesAug2019/checkSNPsJan2020/mapped_reads
name=Adi007
ref=$name.reference.fasta
thr=70
read=$name.fastq

# generate this index before loading Sentieon if you ever want the index to finish
# time ( bwa index $ref )

# module load samtools
# samtools faidx $ref


module load sentieon-genomics/201808.01-opfuvzr

time ( bwa mem -M -K 10000000 -R "@RG\tID:$name \tSM:$name \tPL:ILLUMINA" -t $thr $ref $DIR/$read > $name.sam  ) &> $name.bwa.timelog

time ( sentieon util sort -o $name.sort.bam -t $thr --sam2bam -i $name.sam  ) &> $name.sam2bam.timelog

time ( sentieon driver -t $thr -r $ref -i $name.sort.bam --algo GCBias --summary $name.GC.summary $name.GC.metric --algo MeanQualityByCycle $name.MQ.metric --algo QualDistribution $name.QD.metric --algo InsertSizeMetricAlgo $name.IS.metric --algo AlignmentStat $name.ALN.metric ) &> $name.metric.timelog

time ( sentieon plot metrics -o $name.metric.pdf gc=$name.GC.metric mq=$name.MQ.metric qd=$name.QD.metric isize=$name.IS.metric ) &> $name.metricPlot.timelog

#### these are for removing duplicate reads
time ( sentieon driver -t $thr -i $name.sort.bam --algo LocusCollector --fun score_info $name.score ) &> $name.locusCollect.timelog

time ( sentieon driver -t $thr -i $name.sort.bam --algo Dedup --rmdup --score_info $name.score --metrics $name.dedup.metric $name.dedup.bam ) &> $name.rmDup.timelog
####

time ( sentieon driver -t $thr -r $ref -i $name.dedup.bam --algo Realigner $name.realign.bam ) &> $name.realign.timelog

### same as HaplotypeCaller
time ( sentieon driver -t $thr -r $ref -i $name.realign.bam  --algo Genotyper $name.gVCF --emit_mode gvcf ) &> $name.gvcf.timelog

# time ( sentieon driver -t $thr -r $ref --algo GVCFtyper $name.vcf *.gVCF ) &> $name.vcf.timelog
# module load snphylo
# sed 's/Chr//g' $name.vcf > $name.numOnly.vcf
# snphylo.sh -v $name.numOnly.vcf -c 5 -P $phylo -b -B 10000

