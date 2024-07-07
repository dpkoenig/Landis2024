#Script to plot parental genetic relationships calculated using PLINK
#The input datasets are created from CCII_PARENTS_AND_EXOME.vcf.gz
#Using Plink and vcftools
######################################################################
library(RColorBrewer)
library(data.table)
#Import accessory functions
source("~/bigdata/BARLEY_CCII_FULL/SCRIPTS/R/func.R")

#Set working directory
setwd('~/bigdata/BARLEY_CCII_FULL/DATA/OUTPUT/')

#Import plink PCA analysis
x<-read.table("CALL_EXOME_2/PLINK/plink_file.eigenvec")
y<-read.table("CALL_EXOME_2/PLINK/plink_file.eigenval")

#Get color scheme
parentcols<-getlinecols()

#Generate plot
pdf("FIGURE/PCA_EXOME_PARENTS.pdf")
plot(x[,3:4],pch=21,bg="grey",col="darkgrey",las=2,
     xlab=paste("PC 1 (",round(100*y[1,1]/sum(y[,1]),digits = 1),"%)",sep=""),
     ylab=paste("PC 2 (",round(100*y[2,1]/sum(y[,1]),digits = 1),"%)",sep=""),
     axes=F)
axis(1)
axis(2)
points(x[c(1:28),3:4],
       bg=parentcols[x[1:28,1]],
       pch=21,
       lwd=.5,
       cex=2)
dev.off()

#Read in allele count data in CCII parents in exome dataset
paraf<-as.data.frame(fread("CALL_EXOME_2/par_counts.frq.count",skip=1))
exaf<-as.data.frame(fread("CALL_EXOME_2/exome_counts.frq.count",skip=1))

#Generate allele frequencies and explore only sites variable in the exome dataset
paraf<-paraf[,6]/paraf[,4]
exaf<-exaf[,6]/exaf[,4]
paraf<-paraf[exaf>0]
exaf<-exaf[exaf>0]
length(exaf)#1,316,845

#Calculate fraction of sites segregating
sum(paraf>0,na.rm=T)/length(exaf)#64.8%
sum(paraf[exaf>.1]>0,na.rm=T)/length(exaf[exaf>.1])#96.4%

#Plot fraction segregating by allele frequency
binsex<-cut(exaf,seq(0,1,.05),include.lowest = T)
pdf("FIGURE/FRACTION_EXOME_SEGREGATING.pdf")
barplot(rbind(table(binsex[paraf>0])/table(binsex),
              1-table(binsex[paraf>0])/table(binsex))[1,],
        las=2,
        ylim=c(0,1),
        ylab="Fraction Segregating")
dev.off()

#Compare allele frequencies across the two datasets
parbins<-cut(paraf,
             seq(0,1,.05),
             include.lowest = T)

exbins<-cut(exaf,
             seq(0,1,.05),
             include.lowest = T)
#
parspec<-table(parbins)
exspec<-table(exbins)
#
pdf("FIGURE/Afs_ex_par.pdf")
barplot(rbind(exspec,parspec),beside=T)
dev.off()
#
twodaf<-data.frame(list(ex=exbins,par=parbins))
twod<-table(twodaf)
#
plotit<-melt(twod+1)

pdf("FIGURE/2D_AFS.pdf",height=6, width=8)
ggplot(plotit,aes(ex,par,fill=value)) +
  geom_tile() +
  scale_fill_viridis(trans='log10') +
  xlab("Global allele frequency") +
  ylab("CCII parent allele frequency")
dev.off()

cor(exaf,paraf,use="pairwise",method="spearman")#R=0.81      
Rtestout<-cor.test(exaf,paraf,use="pairwise",method="spearman",exact=FALSE)#rho=0.81      
