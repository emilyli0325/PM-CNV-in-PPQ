#!/bin/bash
echo "mapping ${1}..."

./bwa3/bwa-0.7.15/bwa mem ./Pfal46 ./${1}.R1.fastq.gz ./${1}.R2.fastq.gz -t 12 -M -R "@RG\tID:${1}\tLB:${1}\tPL:ILLUMINA\tPM:HISEQ\tSM:${1}" > SAM/${1}.sam

java -Xmx30g -jar ./picard/picard.jar SortSam \
     INPUT=SAM/${1}.sam \
     OUTPUT=sorted.bam/${1}.sorted.bam \
     SORT_ORDER=coordinate

java -Xmx30g -jar ./picard/picard.jar MarkDuplicates \
     INPUT=sorted.bam/${1}.sorted.bam \
     OUTPUT=dedup.sorted.bam/${1}.dedup.sorted.bam \
     METRICS_FILE=metrics/${1}.metrics.txt

cd dedup.sorted.bam

java -Xmx30g -jar ./picard/picard.jar BuildBamIndex \
     INPUT=${1}.dedup.sorted.bam

cd ..

echo "BQSR $file..."	 
	java -Xmx30g -jar ./GATK/GenomeAnalysisTK.jar\
		 -T BaseRecalibrator\
		 -R ./Pfal32_GATK_index/PlasmoDB-46_Pfalciparum3D7_Genome.fasta\
		 -I dedup.sorted.bam/${1}.dedup.sorted.bam\
		 -knownSites ./Known_sites/3d7_hb3.combined.final.karo.sort.vcf\
		 -knownSites ./Known_sites/7g8_gb4.combined.final.karo.sort.vcf\
		 -knownSites ./Known_sites/hb3_dd2.combined.final.karo.sort.vcf\
		 -o BQSR/${1}.recal.table
		 
	java -Xmx30g -jar ./GATK/GenomeAnalysisTK.jar\
		 -T PrintReads\
		 -R ./Pfal32_GATK_index/PlasmoDB-46_Pfalciparum3D7_Genome.fasta\
		 -I dedup.sorted.bam/${1}.dedup.sorted.bam\
		 -BQSR BQSR/${1}.recal.table\
		 -o recal.bam/${1}.recal.bam 
		 
echo "Variant calling ${1}..."
	java -Xmx30g -jar ./GATK/GenomeAnalysisTK.jar\
         -T HaplotypeCaller\
         -R ./Pfal32_GATK_index/PlasmoDB-46_Pfalciparum3D7_Genome.fasta\
         -I recal.bam/${1}.recal.bam\
         --emitRefConfidence BP_RESOLUTION\
         -o g.vcf.ploidy2/${1}.g.vcf\
         --useNewAFCalculator \
         --sample_ploidy 2 \
         --dontUseSoftClippedBases \
         --num_cpu_threads_per_data_thread 12
        echo "${1} Variant calling done!"

