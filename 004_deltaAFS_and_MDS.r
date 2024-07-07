#Calculate MDS to visualize evolution of the population over time
####################################################################################
setwd('~/bigdata/BARLEY_CCII_FULL/DATA/OUTPUT/')
###########################libraries################################################
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(caTools)
library(viridis)
library(data.table)
library(reshape2)
library(caTools)
library(ggridges)
####################################################################################
####################AFS Analysis####################################################
####################################################################################
source("../../SCRIPTS/R/func.R")
#Data
#Line names directly from the original vcf
#plink used to generate distance matrices
#Chromosome lengths file generated using samtools idxstats on a aligned BAM
ac<-fread("FILTERED/filtered_ac.txt", stringsAsFactors = F)
ac<-as.data.frame(ac,stingsAsFactors=F)
gtnum<-fread("FILTERED/filtered_gtnum.txt", stringsAsFactors = F)
gtnum<-as.data.frame(gtnum,stingsAsFactors=F)
lnames<-read.table("FULL_PARENTAL_VCF/LINE_NAMES.txt",stringsAsFactors = F)[,1]
pruned<-read.table("PLINK/plink_file.prune.in",stringsAsFactors = F)#[,1]
plinkdm<-read.table("PLINK/plink_file_prune.dist")
colnames(plinkdm)<-rownames(plinkdm)<-read.table("PLINK/plink_file_prune.dist.id")[,1]
cclinfo<-read.table("../INPUT/PARENTAL_DATA.txt",stringsAsFactors = F, header=T, na.strings = " -",sep="\t")
chrslen<-read.table("../INPUT/CHRLEN.txt",stringsAsFactors = F)

#getgraphcols
parentcols<-getlinecols()

#get parent and progeny data and convert to AF
ac.cov<-ac[,seq(5,12,2)]+ac[,seq(6,12,2)]
ac.p<-ac[,seq(6,12,2)]/ac.cov
gtnum<-gtnum/2

#combine dataset
ac.i<-do.call(paste, c(ac[,1:2], sep = ":"))
acgt<-cbind(ac.p,gtnum)
colnames(acgt)<-c("F2","F18","F28","F58",lnames)
macgt<-acgt
macgt[macgt[,1]>0.5,]<-1-macgt[macgt[,1]>0.5,]

#################################################################################################################
#Genome wide average Fst between generations 
#
x<-ac.p[,2]
y<-ac.p[,3]
N<-x*((1-y)-(1-x)) + y*((1-x)-(1-y))
D<-x*(1-y) + y*(1-x)
Fst<-N/D
#
pwFst<-c(calcFst(ac.p[,1],ac.p[,2]),
         calcFst(ac.p[,1],ac.p[,3]),
         calcFst(ac.p[,1],ac.p[,4]),
         calcFst(ac.p[,2],ac.p[,3]),
         calcFst(ac.p[,2],ac.p[,4]),
         calcFst(ac.p[,3],ac.p[,4]))
pwFstsq<-matrix(NA,nrow=4,ncol=4)                
pwFstsq[upper.tri(pwFstsq,diag = F)] <-pwFst                
#
fmod<-lm(pwFst~c(18,28,58,33,62,53))
summary(fmod) # p = 0.00383 Adjusted R-squared:  0.8758
fstvgen<-data.frame(Gendelta=c(18,28,58,33,62,53),Fst=pwFst)
#
pdf("FIGURE/Fst_by_time.pdf")
ggplot(fstvgen,aes(Gendelta,Fst)) +
  geom_point() +
  geom_smooth(method='lm', colour="black", linewidth=0.5) + 
  theme_minimal()
dev.off()

#################################################################################################################
#Representation of private alleles over time
priv<-macgt[rowSums(macgt[,-(1:4)]>0,na.rm=T)==1,]
privseg<-rbind(
  colSums(priv[,-(1:4)]>0 & priv[,1]>0,na.rm=T)/colSums(priv[,-(1:4)]>0,na.rm=T),
  colSums(priv[,-(1:4)]>0 & priv[,2]>0,na.rm=T)/colSums(priv[,-(1:4)]>0,na.rm=T),
  colSums(priv[,-(1:4)]>0 & priv[,3]>0,na.rm=T)/colSums(priv[,-(1:4)]>0,na.rm=T),
  colSums(priv[,-(1:4)]>0 & priv[,4]>0,na.rm=T)/colSums(priv[,-(1:4)]>0,na.rm=T))

nprivfix<-colSums(priv[,-(1:4)]>0 & priv[,4]==1,na.rm=T)
sort(nprivfix/sum(priv[,4]==1,na.rm=T))
nminfix<-colSums(macgt[,-(1:4)]>0 & macgt[,4]==1,na.rm=T)
sort(nminfix/sum(macgt[,4]==1,na.rm=T))
#Plot private allele representation over time
pdf("FIGURE/Private_allele_fate.pdf")
plot(0,type="n",ylim=c(0,1),xlim=c(0,60),xlab="Generation",ylab="Fraction of Private Alleles Segregating",axes=F)
abline(v=c(18,28,58))
for(i in 1:nrow(t(privseg))){
  lines(c(2,18,28,58),t(privseg)[i,],col=parentcols[colnames(privseg)[i]],lwd=2)
}
axis(1,las=2)
axis(2,las=2)
dev.off()
#################################################################################################################
#Fix
fix.i<-c(
  (macgt[,2]==1 & 
     macgt[,3]==1 & 
     macgt[,4]==1) |
    (macgt[,3]==1 & 
       macgt[,4]==1) |
    macgt[,4]==1
)
#
macgt_fix<-macgt[fix.i,]          
length(macgt_fix)
colSums(macgt_fix>0,na.rm=T)/colSums(!is.na(macgt_fix))

#################################################################################################################
#Fate of Minor Alleles based on initial allele frequency
F18up<-macgt[,2]==1 & macgt[,3]==1 & macgt[,4]==1 
F28up<-macgt[,3]==1 & macgt[,4]==1 
F58up<-macgt[,4]==1 
F18down<-macgt[,2]==0 & macgt[,3]==0 & macgt[,4]==0 
F28down<-macgt[,3]==0 & macgt[,4]==0 
F58down<-macgt[,4]==0 

#fraction fixed
founding_af<-cut(macgt[,1],seq(0,.5,.05))
plot_frac_fix<-melt(
  data.frame(
    list(
      AF=names(table(founding_af)),
      F18=as.vector(tapply(F18up | F18down,founding_af,sum)/table(founding_af)),
      F28=as.vector(tapply(F28up | F28down,founding_af,sum)/table(founding_af)),
      F58=as.vector(tapply(F58up | F58down,founding_af,sum)/table(founding_af))
      )
    )
)
#
pdf("FIGURE/Fraction_fixedxMAF.pdf")
ggplot(plot_frac_fix,aes(x=AF,y=value,group=variable,fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_viridis(discrete=TRUE,direction = -1) +
  theme_classic()
dev.off()

#fraction fixed up vs down
founding_af<-cut(macgt[,1],seq(0,.5,.05))
plot_fateaf<-melt(
  data.frame(
    list(
      AF=names(tapply(F18up,founding_af,sum)/tapply(F18up | F18down,founding_af,sum)),
      F18=tapply(F18up,founding_af,sum)/tapply(F18up | F18down,founding_af,sum),
      F28=tapply(F28up,founding_af,sum)/tapply(F28up | F28down,founding_af,sum),
      F58=tapply(F58up,founding_af,sum)/tapply(F58up | F58down,founding_af,sum))
    )
)

pdf("FIGURE/Fraction_fixedupxMAF.pdf")
ggplot(plot_fateaf,aes(x=AF,y=value,group=variable,color=variable)) +
  geom_line() +
  scale_color_viridis(discrete=TRUE,direction = -1) +
  theme_classic()
dev.off()

#################################################################################################################
#Fate Of minor alleles carried by parent over time
deltaaf<-macgt-macgt[,1]
large.i<-deltaaf>0
maf2p<-rbind(colSums(macgt[large.i[,1],-(1:4)],na.rm=T)/colSums(!is.na(macgt[,-(1:4)])),
             colSums(macgt[large.i[,2],-(1:4)],na.rm=T)/colSums(!is.na(macgt[,-(1:4)])),
             colSums(macgt[large.i[,3],-(1:4)],na.rm=T)/colSums(!is.na(macgt[,-(1:4)])),
             colSums(macgt[large.i[,4],-(1:4)],na.rm=T)/colSums(!is.na(macgt[,-(1:4)])))

pdf("FIGURE/Fate_of_minor_alleles_by_parents.pdf")
plot(c(0,18,28,58),maf2p[,2],ylim=c(0,max(maf2p)),type="n",xlab="Generation",ylab="Fraction of minor alleles",las=2)
for (i in 1:28){
  lines(c(0,18,28,58),maf2p[,i],col = parentcols[colnames(maf2p)[i]])
}
dev.off()

#################################################################################################################
#MDS Analysis of parents and pooled sequencing AFS
#################################################################################################################
#pruned sites only to generate MDS dataset
prunedpaste<-paste(pruned[,1],pruned[,2],sep=":")
acgt_p<-acgt[ac.i%in%prunedpaste,]

#generate nei distance for parents and progeny pools in the pruned dataset
dm<-matrix(NA,ncol=ncol(acgt_p),nrow=ncol(acgt_p))
for(i in 1:(ncol(dm)-1)){
  for(j in i:ncol(dm)){
    dm[i,j]<-neid(acgt_p[,i],acgt_p[,j])
    dm[j,i]<-dm[i,j]
    print(paste(i,j))
  }
}
colnames(dm)<-rownames(dm)<-colnames(acgt_p)
sort(dm[1,])
sort(dm[2,])
sort(dm[3,])
sort(dm[4,])

#Calculate distance to parents over time
dm2par<-dm[-(1:4),1:4]
dm2par<-dm2par[,1]-dm2par

#
pdf("FIGURE/Distance_to_parent.pdf")
plot(1,type="n",xlim=c(0,60),ylim=c(-.2,.05),las=2,ylab="1-Nei's Genetic Distance",xlab="Generation")
abline(h=0,lty=2)
for(i in 1:nrow(dm2par)){
  lines(c(2,18,28,58),dm2par[i,],col=parentcols[rownames(dm2par)[i]],lwd=2)
}
dev.off()

#Conduct and plot MDS analysis
fit <- cmdscale(as.dist(dm),eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

pdf("FIGURE_R/MDS_with_progeny.pdf")
plot(x, y,
     xlab="Coordinate 1",
     ylab="Coordinate 2",
     main="Metric MDS",
     axes=F,
     type="n",
     xlim=c(-.2,.2),
     ylim=c(-.3,.2))
points(x[rownames(dm)[-(1:4)]],y[rownames(dm)[-(1:4)]],pch=16,co=parentcols[rownames(dm)[-(1:4)]],cex=2)
points(x[1:4],y[1:4],pch=21,bg=brewer.pal(9,"Reds")[c(2,4,6,8)],cex=4,lwd=1.5)
text(x,y,colnames(acgt_p))
axis(1,las=2)
axis(2,las=2)
dev.off()

rm(list=ls())
gc()

