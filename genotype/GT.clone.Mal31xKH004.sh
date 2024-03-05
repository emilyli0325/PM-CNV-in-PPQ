################################################################################################
### MOI
# install packages, with R3.6.3, https://github.com/bahlolab/moimix

library(moimix)
library(SeqArray)
setwd("")

# Converting a VCF to GDS format
seqVCF2GDS("Mal31xKH004.all.clones.01.13.2020.GATK.SNP.filter.vcf", "Mal31xKH004.clone.gds")
my_vcf <-seqOpen("Mal31xKH004.clone.gds")
seqSummary(my_vcf)

###  load data
isolates <- my_vcf
seqSummary(isolates)
sample.id <- seqGetData(isolates, "sample.id")
coords <- getCoordinates(isolates)
head(coords)

# Estimating the BAF matrix
isolate_baf <- bafMatrix(isolates)
class(isolate_baf)
str(isolate_baf)
plot(isolate_baf, "M4.6.C10")

# fws
fws_all <- getFws(isolates)
write.csv(fws_all, "Mal31xKH004.allclone.FWS.csv")

FWS <- read.csv("Mal31xKH004.allclone.FWS.csv", header=TRUE, sep=",", check.names=FALSE)

################################################################################################
################################################################################################
### select single and clean samples
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
   -V Mal31xKH004.all.clones.01.13.2020.GATK.SNP.filter.vcf \
   --sample_file single.list \
   -o Mal31xKH004.clone.singleInf.01.13.2020.vcf 

### '0/1' to non-call
# HetFilter
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
--variant Mal31xKH004.clone.singleInf.01.13.2020.vcf   \
--genotypeFilterExpression "isHet == 1" \
--genotypeFilterName "HetFilter" \
--setFilteredGtToNocall \
--out Mal31xKH004.clone.singleInf.rmHet.01.13.2020.vcf


java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V Mal31xKH004.clone.singleInf.rmHet.01.13.2020.vcf \
-F CHROM -F POS -F REF -F ALT \
-GF GT -GF AD -GF DP -GF GQ -GF PL \
-o Mal31xKH004.clone.singleInf.rmHet.01.13.2020.table 

################################################################################################
################################################################################################
## hmmIBD

## data clean
Input <- read.delim("Mal31xKH004.clone.singleInf.rmHet.01.13.2020.table", sep="\t",header=T)
> dim(Input)
[1] 13380  3949

GT <- Input[,seq(5,3949,5)]
> dim(GT)
[1] 13380   789


GT.A <- function(X){unlist(strsplit(as.character(X),"\\/"))[1]}

GT.Aallele <- matrix(ncol=789,nrow=13380)
for(j in 1:789){
GT.Aallele[,j] <- sapply(GT[,j], GT.A, simplify="array")
}

GT.Aallele[GT.Aallele== "."]<-NA
colnames(GT.Aallele)<- colnames(GT)
> dim(GT.Aallele)
[1] 13380   789

GT.recode<-NULL
Mal31 <- as.character(GT.Aallele[,789])
KH004 <- as.character(GT.Aallele[,784])

for(i in c(1:789)){
    X<-as.character(GT.Aallele[,i])
    Z<-NULL; 
	Z[which(X==Mal31 & is.na(X)==F)]<- 0; 
	Z[which(X==KH004 & is.na(X)==F)]<- 1; 
	Z[which(is.na(X)==T)]<- "-1";
    GT.recode <- cbind(GT.recode,Z)
}


GT.recode <- cbind(Input[,1:2],GT.recode)

colnames(GT.recode) <- c("chrom", "pos",gsub("\\.GT|\\.\\.GT","",colnames(GT)[1:789]))

GT.recode.cleanChr <- GT.recode
GT.recode.cleanChr$chrom <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",GT.recode$chrom)

write.table(GT.recode.cleanChr,"Mal31xKH004.clone.singleInf.rmHet.01.13.2020.txt",row.names=F,col.names=T,sep="\t", quote=FALSE)

#### hmmIBD
/master/xli/software/hmmIBD/hmmIBD-master/hmmIBD -i Mal31xKH004.clone.singleInf.rmHet.01.13.2020.txt -o Mal31xKH004.787.singleclone.hmmIBD

#### unique recombianant
/master/xli/software/hmmIBD/hmmIBD-master/hmmIBD -i Mal31xKH004.clone.singleInf.rmHet.01.13.2020.txt -o Mal31xKH004.unique.hmmIBD -b Non.uniqueID.txt
