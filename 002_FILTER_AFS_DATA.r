#Script to filter Allele frequency and parental genotype calls       #
######################################################################
setwd('~/bigdata/BARLEY_CCII_FULL/DATA/OUTPUT/')
library(data.table)

#Read In Data
#gt file generated from FILTERED_PARENTAL_SNPS.vcf to contain only numeric (0,1,2 homref,het,homalt and NA for missing data)
gtnum<-as.data.frame(fread("FULL_PARENTAL_VCF/GT_ONLY_NUMERIC.txt",stringsAsFactors = F))
ac<-as.data.frame(fread("CALL_POOL_AFS/FINAL_AFS.txt",stringsAsFactors = F))

#Final site filter
filt.i<-rowSums(ac[,5:12])<300 & rowSums((ac[,seq(5,12,2)]+ac[,seq(6,12,2)])>=10)==4 & ac[,5]>0 & ac[,6]>0
gtnum<-gtnum[filt.i,]
ac<-ac[filt.i,]

#write downsampled dataset
dir.create("FILTERED")
write.table(ac,"FILTERED/filtered_ac.txt",quote=F,row.names=F,col.names=F)
write.table(gtnum,"FILTERED/filtered_gtnum.txt",quote=F,row.names=F,col.names=F)
rm(list = ls())
gc()

