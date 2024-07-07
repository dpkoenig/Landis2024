#Import and plot the selection scan results
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
####################################################################################
####################################################################################
source("../../SCRIPTS/R/func.R")
#read in data
ac<-fread("FILTERED/filtered_ac.txt",stringsAsFactors = F)
ac<-as.data.frame(ac)

#calculate allele frequency change over time
ac.cov<-ac[,seq(5,12,2)]+ac[,seq(6,12,2)]
ac.p<-ac[,seq(6,12,2)]/ac.cov

#################################################################################################################
#Windowed change in af
maf.ac.p<-ac.p
maf.ac.p[maf.ac.p[,1]>.5,]<-1-maf.ac.p[maf.ac.p[,1]>.5,]
delta.maf<-maf.ac.p-maf.ac.p[,1]
win.delta.maf<-caTools::runmean(as.matrix(delta.maf),1000,align = "left")
win.delta.maf<-win.delta.maf[seq(1,nrow(win.delta.maf),500),]
del.af<-abs(ac.p-ac.p[,1])
win.delta.af<-caTools::runmean(as.matrix(del.af),1000,align = "left")
win.delta.af<-win.delta.af[seq(1,nrow(win.delta.af),500),]

#Windowed expected heterozygozity
He<-2*ac.p*(1-ac.p)
Hewins<-caTools::runmean(as.matrix(He),1000,align = "left")
Hewins<-Hewins[seq(1,nrow(Hewins),500),]
startpos<-runmin(ac[,2],1000,align = "left")
startpos<-startpos[seq(1,length(startpos),500)]
endpos<-runmax(ac[,2],1000,align = "left")
endpos<-endpos[seq(1,length(endpos),500)]
chr<-ac[seq(1,nrow(ac),500),1]

#
colnames(Hewins)<-c("F0",
                    "F18",
                    "F28",
                    "F58")
colnames(win.delta.af)<-c("F0",
                          "F18",
                          "F28",
                          "F58")
colnames(win.delta.maf)<-c("F0",
                           "F18",
                           "F28",
                           "F58")

#####################################################################################
#Density plot of He and del.af over time
#####################################################################################
Hewins.m<-melt(Hewins)
win.delta.af.m<-melt(win.delta.af[,-1])
win.delta.maf.m<-melt(win.delta.maf[,-1])
colnames(Hewins.m)<-c("Var1","Generation","He")
colnames(win.delta.af.m)<-c("Var1","Generation","daf")
colnames(win.delta.maf.m)<-c("Var1","Generation","dmaf")

#
pdf("FIGURE/Change_in_AF_He_windowed.pdf")
ggplot(Hewins.m,aes(x=He,color=Generation)) +
  geom_density(lwd=1.5,bw=.01) +
  scale_color_viridis(discrete=TRUE,direction = -1) +
  theme_classic()

ggplot(win.delta.af.m,aes(x=daf,color=Generation)) +
  geom_density(lwd=1.5,bw=.01) +
  scale_color_viridis(end=.75,discrete=TRUE,direction = -1) +
  theme_classic()

ggplot(win.delta.maf.m,aes(x=dmaf,color=Generation)) +
  geom_density(lwd=1.5,bw=.02) +
  scale_color_viridis(end=.75,discrete=TRUE,direction = -1) +
  theme_classic()
dev.off()

#####################################################################################
#Calculating and plotting AFS
#####################################################################################
#####################################################################################
#1D
maf<-rbind(table(cut(maf.ac.p[,1],seq(0,1,.05),include.lowest = T)),
           table(cut(maf.ac.p[,2],seq(0,1,.05),include.lowest = T)),
           table(cut(maf.ac.p[,3],seq(0,1,.05),include.lowest = T)),
           table(cut(maf.ac.p[,4],seq(0,1,.05),include.lowest = T)))
rownames(maf)<-c("P",
                 "F18",
                 "F28",
                 "F58")
afsspec<-melt(maf)
colnames(afsspec)<-c("Generation","AF","Count")

#Plot
pdf("FIGURE/AFS_1D_overtime.pdf")
ggplot(afsspec,aes(x=AF,y=Count,fill=Generation)) +
  geom_bar(stat= "identity",position='dodge') +
  scale_fill_viridis(discrete=T,direction=-1) +
  theme_minimal()
dev.off()

#####################################################################################
#2d
for2d<-list()
for (i in 1:ncol(maf.ac.p)){
  for2d[[i]]<-cut(maf.ac.p[,i],seq(0,1,.05),include.lowest=T)
}
for2d<-as.data.frame(for2d)
colnames(for2d)<-c("P",
                   "F18",
                   "F28",
                   "F58")

#2D AFS Plot
pdf("FIGURE/2d_AFS.pdf")
for (i in 1:3){
  for (j in (i+1):4){
    the2d<-table(for2d[,3],for2d[,4])
    the2d_A_B<-melt(the2d)
    ggplot(the2d_A_B,aes(x=Var1,y=Var2,fill=value)) +
      geom_tile() +
      scale_fill_viridis(trans = "log")
  }
}
dev.off()

####################################################################################
#Selection Scan
#####################################################################################

#####################################################################################
#Create data frame
finalpos<-data.frame(list(chr=chr,
                          startpos=startpos,
                          endpos=endpos,
                          He0=Hewins[,1],
                          He18=Hewins[,2],
                          He28=Hewins[,3],
                          He58=Hewins[,4],
                          daf=win.delta.af)
)
finalpos<-finalpos[!finalpos$chr%in%"chrUn",]
finalpos<-finalpos[finalpos$endpos-finalpos$startpos<10000000,]

#
nearfix<-c(sum(width(reduce(GRanges(finalpos[finalpos[,4]<.01,1],IRanges(finalpos[finalpos[,4]<.01,2],finalpos[finalpos[,4]<.01,3])))))/sum(width(reduce(GRanges(finalpos[,1],IRanges(finalpos[,2],finalpos[,3]))))),
  sum(width(reduce(GRanges(finalpos[finalpos[,5]<.01,1],IRanges(finalpos[finalpos[,5]<.01,2],finalpos[finalpos[,5]<.01,3])))))/sum(width(reduce(GRanges(finalpos[,1],IRanges(finalpos[,2],finalpos[,3]))))),
  sum(width(reduce(GRanges(finalpos[finalpos[,6]<.01,1],IRanges(finalpos[finalpos[,6]<.01,2],finalpos[finalpos[,6]<.01,3])))))/sum(width(reduce(GRanges(finalpos[,1],IRanges(finalpos[,2],finalpos[,3]))))),
  sum(width(reduce(GRanges(finalpos[finalpos[,7]<.01,1],IRanges(finalpos[finalpos[,7]<.01,2],finalpos[finalpos[,7]<.01,3])))))/sum(width(reduce(GRanges(finalpos[,1],IRanges(finalpos[,2],finalpos[,3]))))))
#0.000000000 0.007001123 0.008139874 0.298007401

#Read in simulation data for each chromosome region
onehunN<-read.table("../INPUT/combined_100N_100000S.txt",header=F,stringsAsFactors = F)

#Order and remove unmapped contigs
onehunN<-onehunN[order(onehunN[,1],onehunN[,2]),]
onehunN<-onehunN[!onehunN[,1]%in%"chrUn",]
onehunN<-onehunN[onehunN[,2]<onehunN[,4],]

#Create final combined dataset
finalpos<-cbind(finalpos,dAFsim=onehunN[,6],Hesim=onehunN[,7])

#####################################################################################
#calculate pvalues for scan
#####################################################################################
pvalHe<-finalpos$Hesim/100001
pvaldAF<-(100002-finalpos$dAFsim)/100001
chiall<--2*(log(2*pvaldAF)+log(pvalHe))
bonthresh<--log(.05/nrow(finalpos),10)
pval<-1-(pchisq(chiall,4))
sig.i<-pvalHe<.05 & pvaldAF<.05 & p.adjust(pval,method="BH")<.05# & finalpos$He18<.05 & finalpos$He28<.1 & finalpos$He58<.1

# Collapse Get significant regions
sigforpub<-as.data.frame(reduce(IRanges((1:nrow(finalpos))[sig.i],(1:nrow(finalpos))[sig.i])+1)-1)
sigforpub<-sigforpub[sigforpub[,3]>1,]
sigreg<-as.data.frame(findOverlaps(IRanges(1:nrow(finalpos),1:nrow(finalpos)),IRanges(sigforpub[,1],sigforpub[,2])))[,1]
forpub<-reduce(GRanges(finalpos[sigreg,1],IRanges(finalpos[sigreg,2],finalpos[sigreg,3])))

# Generate final scan matrix
selscanfinal<-data.frame(finalpos,pvalHe,pvaldAF,pval,sig.i)
sum(as.data.frame(forpub)[,4])/sum(tapply(selscanfinal$endpos,selscanfinal$chr,max))

# Write out full matrix of selection scan data and the significant regions
write.table(selscanfinal,"TABLE/final_sel_scan.txt",quote=F)
write.table(as.data.frame(forpub),"TABLE/final_selected_regions.txt",quote=F)

###############################################################################################
# Plot selection scan
chrcol<-rep(brewer.pal(8,"Set2")[8],nrow(finalpos))
chrcol[finalpos[,1]%in%c("chr2H","chr4H","chr6H")]<-brewer.pal(8,"Pastel2")[8]
chrcol2<-rep(brewer.pal(8,"Set2")[3],nrow(finalpos))
#chrcol2[finalpos[,1]%in%c("chr2H","chr4H","chr6H")]<-brewer.pal(8,"Set2")[3]

#Plot Selection Scan
pvalsig<--log(pval,10)
pvalsig[-sigreg]<-NA
pdf("FIGURE/sel_scan.pdf",height=4,width=12)
plot(-log(pval,10),pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,8))
points(pvalsig,cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
axis(2,las=2)
dev.off()

newpos<-rowMeans(finalpos[,2:3])
for(i in 2:7){
  thischr<-unique(finalpos[,1])[i]
  lastchr<-unique(finalpos[,1])[i-1]
  newpos[finalpos[,1]==thischr]<-newpos[finalpos[,1]==thischr]+100000000+max(newpos[finalpos[,1]==lastchr])
}

#
pdf("FIGURE/sel_scan_w_pos.pdf",height=4,width=8)
plot(newpos,-log(pval,10),pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,8))
points(newpos,pvalsig,cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
axis(2,las=2)
axis(1,las=2,at=tapply(newpos,finalpos[,1],min),labels = rep(0,7),lwd=0,lwd.ticks = 1)
axis(1,las=2,at=tapply(newpos,finalpos[,1],max),labels = round(tapply(rowMeans(finalpos[,2:3]),finalpos[,1],max)/1e6,2),lwd=0,lwd.ticks = 1)
dev.off()

###################################################################################################
pdf("FIGURE/He_and_daf_scan_w_pos.pdf",height=10,width=7.5)
###################################################################################################
###################################################################################################
Hesig<-finalpos$He18
Hesig[-sigreg]<-NA
par(mfrow=c(7,1),mai=c(0.1,0.1,0.1,0.1),oma=c(4,4,4,4))
for (i in unique(finalpos[,1])){
  plot(rowMeans(finalpos[finalpos$chr==i,2:3]),finalpos$He18[finalpos$chr==i],pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,0.5),xlim=c(0,max(rowMeans(finalpos[,2:3]))))
  points(rowMeans(finalpos[finalpos$chr==i,2:3]),Hesig[finalpos$chr==i],cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
  axis(2,las=2)
}
axis(1,las=2)

###################################################################################################
###################################################################################################
Hesig<-finalpos$He28
Hesig[-sigreg]<-NA
par(mfrow=c(7,1),mai=c(0.1,0.1,0.1,0.1),oma=c(4,4,4,4))
for (i in unique(finalpos[,1])){
  plot(rowMeans(finalpos[finalpos$chr==i,2:3]),finalpos$He28[finalpos$chr==i],pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,0.5),xlim=c(0,max(rowMeans(finalpos[,2:3]))))
  points(rowMeans(finalpos[finalpos$chr==i,2:3]),Hesig[finalpos$chr==i],cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
  axis(2,las=2)
}
axis(1,las=2)

###################################################################################################
###################################################################################################
Hesig<-finalpos$He58
Hesig[-sigreg]<-NA
par(mfrow=c(7,1),mai=c(0.1,0.1,0.1,0.1),oma=c(4,4,4,4))
for (i in unique(finalpos[,1])){
  plot(rowMeans(finalpos[finalpos$chr==i,2:3]),finalpos$He58[finalpos$chr==i],pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,0.5),xlim=c(0,max(rowMeans(finalpos[,2:3]))))
  points(rowMeans(finalpos[finalpos$chr==i,2:3]),Hesig[finalpos$chr==i],cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
  axis(2,las=2)
}
axis(1,las=2)

###################################################################################################
###################################################################################################
dafsig<-finalpos$daf.F18
dafsig[-sigreg]<-NA
par(mfrow=c(7,1),mai=c(0.1,0.1,0.1,0.1),oma=c(4,4,4,4))
for (i in unique(finalpos[,1])){
  plot(rowMeans(finalpos[finalpos$chr==i,2:3]),finalpos$daf.F18[finalpos$chr==i],pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,0.7),xlim=c(0,max(rowMeans(finalpos[,2:3]))))
  points(rowMeans(finalpos[finalpos$chr==i,2:3]),dafsig[finalpos$chr==i],cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
  axis(2,las=2)
}
axis(1,las=2)

###################################################################################################
###################################################################################################
dafsig<-finalpos$daf.F28
dafsig[-sigreg]<-NA
par(mfrow=c(7,1),mai=c(0.1,0.1,0.1,0.1),oma=c(4,4,4,4))
for (i in unique(finalpos[,1])){
  plot(rowMeans(finalpos[finalpos$chr==i,2:3]),finalpos$daf.F28[finalpos$chr==i],pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,0.7),xlim=c(0,max(rowMeans(finalpos[,2:3]))))
  points(rowMeans(finalpos[finalpos$chr==i,2:3]),dafsig[finalpos$chr==i],cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
  axis(2,las=2)
}
axis(1,las=2)

###################################################################################################
###################################################################################################
dafsig<-finalpos$daf.F58
dafsig[-sigreg]<-NA
par(mfrow=c(7,1),mai=c(0.1,0.1,0.1,0.1),oma=c(4,4,4,4))
for (i in unique(finalpos[,1])){
  plot(rowMeans(finalpos[finalpos$chr==i,2:3]),finalpos$daf.F58[finalpos$chr==i],pch=16,cex=.5,col=chrcol,lwd=.1,axes=F,xlab="",ylim=c(0,0.7),xlim=c(0,max(rowMeans(finalpos[,2:3]))))
  points(rowMeans(finalpos[finalpos$chr==i,2:3]),dafsig[finalpos$chr==i],cex=.7,pch=21,bg=chrcol2,col=rgb(0,0,0,0.9),lwd=.1)
  axis(2,las=2)
}
axis(1,las=2)
dev.off()

#
bychrAll<-split(-log(pval,10),finalpos[,1])
bychrSig<-split(pvalsig,finalpos[,1])
bychrpos<-split(rowMeans(finalpos[,2:3]),finalpos[,1])

i<-2
plot(bychrpos[[i]],bychrAll[[i]],pch=16,cex=.5,col="grey",lwd=.1,axes=F,xlab="",ylim=c(0,8))
points(bychrpos[[i]],bychrSig[[i]],cex=.7,pch=21,bg=chrcol2,col="black",lwd=.1)
axis(1,las=2)
axis(2,las=2)

###############################################################################################
#Compare significant vs background He
forplot<-melt(selscanfinal[,c(4:7,17)])

#Plot
pdf("FIGURE/significant_Reg_He.pdf")
ggplot(forplot,aes(variable,value)) +
  geom_boxplot(aes(color=sig.i),outlier.size = .5) +
  scale_color_viridis(discrete = T,end=.3) +
  theme_classic()
dev.off()

t.test(selscanfinal$He0~sig.i)$p.value # 1.827377e-187
t.test(selscanfinal$He18~sig.i)$p.value # 0
t.test(selscanfinal$He28~sig.i)$p.value # 8.85518e-95
t.test(selscanfinal$He58~sig.i)$p.value # 2.331397e-319

###############################################################################################
#Compare significant vs background daf
forplot<-melt(selscanfinal[,c(9:11,17)])
pdf("FIGURE/significant_Reg_Daf.pdf")
ggplot(forplot,aes(variable,value)) +
  geom_boxplot(aes(color=sig.i),outlier.size = .5) +
  scale_color_viridis(discrete = T,end=.3) +
  theme_classic()
dev.off()

t.test(selscanfinal$daf.F18~sig.i)$p.value # 1.778636e-322
t.test(selscanfinal$daf.F28~sig.i)$p.value # 9.331963e-293
t.test(selscanfinal$daf.F58~sig.i)$p.value # 2.213878e-243

#Plot whole genome stats
pdf("FIGURE/scans_other_stats.pdf",height=3,width=5)
par(mfrow=c(2,1),mai=c(0.5,1,0.5,0.5))
#
plot(selscanfinal$He0 - selscanfinal$He18,
     col=chrcol,pch=16,ylim=c(0,.5),cex=.7,axes=F,xlab="",ylab="dHe",type="n")
abline(h=seq(0,.5,.05),lty=2,lwd=.5)
points(selscanfinal$He0 - selscanfinal$He18,col=chrcol,pch=16,cex=.5)
tmp<-selscanfinal$He0 - selscanfinal$He18
tmp[!selscanfinal$sig.i]<-NA
points(tmp,pch=16,col="red",cex=.5)
axis(2,las=2)

#
plot(selscanfinal$daf.F18,col=chrcol,pch=16,cex=.7,ylim=c(0,1),axes=F,xlab="",ylab="dAF",type="n")
abline(h=seq(0,1,.1),lty=2,lwd=.5)
points(selscanfinal$daf.F18,col=chrcol,pch=16,cex=.5)
tmp<-selscanfinal$daf.F18
tmp[!selscanfinal$sig.i]<-NA
points(tmp,pch=16,col="blue",cex=.5)
axis(2,las=2)
dev.off()

#breakout chromosomes
pdf("FIGURE/selscanxchr.pdf",height=11,width=8.5)
par(mfrow=c(7,1),mai=c(0.2,0.5,0.2,0.5),oma=c(2,2,2,2))
for (i in unique(selscanfinal$chr)){
  thischrpos<-rowMeans(selscanfinal[selscanfinal$chr==i,2:3])
  thischrpval<-pval[selscanfinal$chr==i]
  thischrpvalsig<-pvalsig[selscanfinal$chr==i]
  plot(thischrpos,-log(thischrpval,10),pch=16,cex=.4,col="darkgrey",axes=F,xlab="",ylim=c(0,8),xlim=c(0,max(selscanfinal[,3])))
  points(thischrpos,thischrpvalsig,cex=.7,pch=21,bg="grey35",lwd=.2)
  axis(2,las=2)
}
axis(1,las=2)
dev.off()
