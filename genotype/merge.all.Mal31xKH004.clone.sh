## g.vcf cohort
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs -R /./PlasmoDB-46_Pfalciparum3D7_Genome.fasta \
-L ./Known_sites/Core_Genome.intervals \
-V ./Mal31.di.g.vcf \
-V ./KH004-2-019-H9.di.g.vcf \
...
--useNewAFCalculator \
--sample_ploidy 2 \
-nt 20 \
-o ./Mal31xKH004.all.clones.vcf

## MT
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V Mal31xKH004.all.clones.vcf \
-F CHROM -F POS -F REF -F ALT \
-GF GT -GF AD \
-o Mal31xKH004.all.clones.01.13.2020.GATK.di.MT.table 

# Seperate SNP and Indel
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V Mal31xKH004.all.clones.vcf \
-selectType SNP -o Mal31xKH004.all.clones.01.13.2020.GATK.SNP.di.vcf


# find SNPs called by both parents
java -Xmx50g -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
   -V Mal31xKH004.all.clones.01.13.2020.GATK.SNP.di.vcf \
   --concordance ./P.genotype_MAL31xKH004/Mal31_KH004.core.SNP.filter.vcf \
   --restrictAllelesTo BIALLELIC -nt 10 \
   -o Mal31xKH004.all.clones.01.13.2020.GATK.SNP.filter.vcf
