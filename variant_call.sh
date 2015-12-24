#!/bin/sh

#run this progam in the analysis folder
PREFIX=$1
DECOY=/data/BIMAS_LMS/Exome/hs37d5/genome.fa.gz
bamDir=../bam2
fastq_dir=../fastq_dir
mdup_metrics=../mdup_metrics
log_dir=../log_dir
tmp_dir=../tmp_dir


module load samtools/0.1.19
module load bwa/0.7.8
module load biobambam
#mkdir mdup_metrics
#mkdir fastq_dir
#mkdir log_dir
bamtofastq exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename=$bamDir/${PREFIX}.bam T=tmp_fq_$PREFIX S=$fastq_dir/${PREFIX}_single.fq O=$fastq_dir/${PREFIX}_unmatched_r1.fq O2=$fastq_dir/${PREFIX}_unmatched_r2.fq collate=1 tryoq=1 > $fastq_dir/${PREFIX}.fq 2>$log_dir/log_bam2fq_$PREFIX
samtools view -H $bamDir/${PREFIX}.bam|grep '^@RG' |sed 's/\t/\\t/g' > ${PREFIX}.RG

bwa mem -t 8 -p -T 0 -R "`cat ${PREFIX}.RG`" $DECOY $fastq_dir/${PREFIX}.fq - |bamsort inputformat=sam level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=$DECOY tmpfile=tmp_bwa_${PREFIX} O=${PREFIX}.sort.bam
rm $fastq_dir/${PREFIX}.fq
rm ${PREFIX}.RG

bammarkduplicates O=${PREFIX}.mdup.bam tmpfile=tmp_mdup_${PREFIX} markthreads=8 rewritebam=1 rewritebamlevel=1 index=1 md5=1 M=$mdup_metrics/mdup_metrics_${PREFIX} I=${PREFIX}.sort.bam 
rm ${PREFIX}.sort.bam
touch ${PREFIX}.pcap.finish

gatk_RES=/fdb/GATK_resource_bundle/b37-2.8
module load GATK/3.3-0
INTERVAL=/data/BIMAS_LMS/Exome/SureSelect_Exon_V5.interval_list
java -Xmx16g -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $gatk_RES/human_g1k_v37_decoy.fasta -I ${PREFIX}.mdup.bam -o ${PREFIX}.realign.intervals -known $gatk_RES/Mills_and_1000G_gold_standard.indels.b37.vcf -known $gatk_RES/1000G_phase1.indels.b37.vcf -L ${INTERVAL} --interval_padding 100 -nt 8
java -Xmx4g -Djava.io.tmpdir=$tmp_dir -jar $GATK_HOME/GenomeAnalysisTK.jar -R $gatk_RES/human_g1k_v37_decoy.fasta -T IndelRealigner -I ${PREFIX}.mdup.bam -targetIntervals ${PREFIX}.realign.intervals -o ${PREFIX}.mdup.realigned.bam  -known $gatk_RES/Mills_and_1000G_gold_standard.indels.b37.vcf -known $gatk_RES/1000G_phase1.indels.b37.vcf --maxReadsInMemory 400000 -maxReads 5000

samtools sort -m 4000000000 ${PREFIX}.mdup.realigned.bam ${PREFIX}.mdup.realigned.sorted
samtools index ${PREFIX}.mdup.realigned.sorted.bam
rm ${PREFIX}.mdup.realigned.bam
rm ${PREFIX}.mdup.realigned.bai
touch ${PREFIX}.mdup.realign.sort.finish

java -Xmx8g -jar $GATK_HOME/GenomeAnalysisTK.jar -l INFO -R $gatk_RES/human_g1k_v37_decoy.fasta -knownSites $gatk_RES/dbsnp_138.b37.vcf -knownSites $gatk_RES/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites $gatk_RES/1000G_phase1.indels.b37.vcf -I ${PREFIX}.mdup.realigned.sorted.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${PREFIX}.recal_data.table -nct 4

java -Xmx8g -jar $GATK_HOME/GenomeAnalysisTK.jar -T PrintReads -R $gatk_RES/human_g1k_v37_decoy.fasta -I ${PREFIX}.mdup.realigned.sorted.bam -BQSR ${PREFIX}.recal_data.table -o ${PREFIX}.mdup.realigned.sorted.recal.bam -nct 4
rm ${PREFIX}.recal_data.table 
rm ${PREFIX}.mdup.realigned.sorted.bam
rm ${PREFIX}.mdup.realigned.sorted.bam.bai
touch ${PREFIX}.mdup.realigned.sorted.recal.bam.finish

:'
#module load bedtools
#BED=/data/BIMAS_LMS/Exome/SureSelect_Exon_V5_+-100.bed
#intersectBed -abam ${PREFIX}.recal.bam -b $BED -wa|samtools view -h -b -q 5 -F 1024 - > ${PREFIX}.exon_+-100.bam
#samtools mpileup -BQ0 -f human_g1k_v37_decoy.fasta ${PREFIX}.realigned.sorted.recal.bam > ${PREFIX}.mpileup
#java -jar /usr/local/apps/varscan/2.3.6/VarScan.v2.3.6.jar  pileup2snp ${PREFIX}.mpileup --p-value 0.98 --min-coverage 2 --min-var-freq 0.1 > ${PREFIX}.varscan.txt
#rm ${PREFIX}.mpileup
#touch ${PREFIX}.varscan.finish
'

