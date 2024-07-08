setwd('~/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/')
##################################################################################################################
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(fields)
library(caTools)
library(GenomicRanges)
library(viridis)
##################################################################################################################
getconcen<-function(x){
  x<-na.omit(x)
  if(is.null(x)){
    return(NA)
  }else{
    outit<-c(sum(x==0,na.rm=T),
             sum(x==1,na.rm=T),
             sum(x==2,na.rm=T))
    return(c(0,1,2)[order(outit,decreasing = T)][1])
  }
}

processg<-function(x){
  if(nrow(x)>1){
    afwi<-colSums(x,na.rm=T)/(2*(colSums(!is.na(x))))
    consen<-rep(NA,length(afwi))
    consen[!is.na(afwi) & afwi<=.5]<-0
    consen[!is.na(afwi) & afwi>.5]<-2
    return(consen)
  }else{
    return(as.vector(x))
  }
}
##################################################################################################################
#Read in datasets
gend<-fread("STITCH/FULL_FILTER.gt",stringsAsFactors = F)
gend<-data.frame(gend,stringsAsFactors=F)[,-1]
colnames(gend)<-read.table("STITCH/FULL_FILTER.gt.names")[,1]
gend.sites<-read.table("STITCH/FULL_FILTER.gt.pos")
happ<-data.frame(fread("../INPUT/plink_geno.fam",stringsAsFactors = F))

##################################################################################################################
##################################################################################################################
#Get distance as fraction
dm<-read.table("PLINK/plink_string.dist")/(2*nrow(gend))
colnames(dm)<-rownames(dm)<-read.table("PLINK/plink_string.dist.id")[,1]
##################################################################################################################
##################################################################################################################
gen.i<-c(rep(0,28),sapply(strsplit(colnames(dm),"_")[-(1:28)],"[",1))
gen.i[gen.i==1]<-18
gen.i[gen.i==2]<-28
gen.i[gen.i==3]<-50
gen.i[gen.i==7]<-58
dm.gen<-list(F0=dm[gen.i==0,gen.i==0],
             F18=dm[gen.i==18,gen.i==18],
             F28=dm[gen.i==28,gen.i==28],
             F50=dm[gen.i==50,gen.i==50],
             F58=dm[gen.i==58,gen.i==58])
forplot<-melt(lapply(dm.gen,function(x){x[upper.tri(x,diag=F)]}))
#
pdf("FIGURES_FILES_TABLES/distance_plot.pdf",height=5,width=3)
ggplot(forplot,aes(value,color=L1)) +
  geom_density(show.legend = F,bw=.025, lwd=1) +
  scale_color_viridis(discrete=TRUE,direction = -1) +
  #  scale_fill_viridis(discrete=TRUE,direction = -1,alpha=.2) +
  theme_classic()
dev.off()

dm2<-dm[-(1:28),-(1:28)]
dm2[lower.tri(dm2,diag = T)]<-NA
fam <-  sapply(strsplit(colnames(dm2),"_"),
               "[",
               1)
dm3<-data.frame(matrix(NA,nrow=nrow(dm2),ncol=ncol(dm2)))
colnames(dm3)<-rownames(dm3)<-rownames(dm2)
for (i in as.character(c(1,2,3,7))){
  dm3[fam==i,fam==i]<-dm2[fam==i,fam==i]
}
#
pdf("FIGURES_FILES_TABLES/distance_matrix.pdf",height=6,width=6)
image.plot(as.matrix(dm3),col=viridis(20),axes=F)
dev.off()
##################################################################################################################
#Cluster Clonal lineages
set.seed(41298371)
dm.prog<-dm
theclust<-hclust(as.dist(dm.prog),method = "average" )
theclust<-cutree(theclust,h=.001)
theclust<-theclust[-(1:28)]
write.table(theclust,"FIGURES_FILES_TABLES/hap_assign.txt")
#
length(table(theclust))#261
sum(table(theclust)==1)#177
sum(table(theclust)==1)/length(table(theclust))#0.6781609

#
forplot<-as.data.frame(table(theclust))
pdf("FIGURES_FILES_TABLES/Clone_Abundance.pdf",height=2,width=5)
ggplot(forplot,aes(x=Freq)) +
  geom_bar() +
  scale_y_continuous(trans='sqrt') +
  theme_classic()
dev.off()

#Get Clone Frequency over time
fam<-sapply(strsplit(names(theclust),"_"),"[",1)
clonexg<-sapply(split(factor(theclust),fam),table)
clonexgfq<-t(t(clonexg)/colSums(clonexg))
clonexgfq<-cbind(1/10000,clonexgfq)

#Number of lineages sampled x generation
colSums(clonexgfq>0)# 93  93  38  55 
apply(clonexgfq,2,max)
clonexgfq["111",]# 0.004587156 0.118181818 0.681818182 0.586363636 
apply(clonexgfq,2,max)#  0.1192661 0.1181818 0.6818182 0.5863636 

chisq.test(cbind(clonexg["111",],
                 colSums(clonexg))[3:4,])#0.3642 p.val
chisq.test(cbind(clonexg["111",],
                 colSums(clonexg)))$p.value#3.781748e-36 p.val
hapdiv<-(1-colSums(clonexgfq[,-1]^2))*(colSums(clonexg)/(colSums(clonexg)-1))# 0.9678688 0.9726858 0.5296389 0.6503944 

#Construct Muller plot
clonexgfq<-clonexgfq[order(rowSums(clonexgfq),decreasing = T),]
getcommon<-clonexgfq[1:10,]
clonexgcumfq<-apply(getcommon,2,cumsum)
clonexgcumfq<-rbind(c(0,0,0,0,0),clonexgcumfq,c(1,1,1,1,1))
colplot<-brewer.pal(nrow(clonexgcumfq)-2,"Spectral")
colplot<-c(rev(colplot),"lightgrey")

#Output Plot
pdf("FIGURES_FILES_TABLES/MULLER.pdf")
plot(0,type="n",xlim=c(2,58),ylim=c(0,1),xaxs="i",yaxs="i",las=2,xlab="Generation",ylab="Relative Frequency")
for (i in 1:(nrow(clonexgcumfq)-1)){
  polygon(c(c(2,18,28,50,58),rev(c(2,18,28,50,58))),c(clonexgcumfq[i,],rev(clonexgcumfq[i+1,])),col=colplot[i],border = "white",)
}
dev.off()

#Get consensus for each clonal lineage
gethaps<-split.data.frame(t(gend[,-(1:28)]),theclust)

#
sum_haps<-lapply(gethaps,processg)
finalhaps<-do.call("cbind",sum_haps)
finalhaps<-apply(finalhaps,2,as.numeric)
write.table(finalhaps,"FIGURES_FILES_TABLES/final_hap_all.txt")
#lin1 is cluster 111 
##################################################################################################################
##################################################################################################################
##################################################################################################################
#Segmental analysis of genetic relationships of lineages to the progeny
finalhaps<-read.table("FIGURES_FILES_TABLES/final_hap_all.txt")
#Generate  object for windowed analysis
pos.i<-cut(gend.sites[,2],seq(1,max(gend.sites[,2]),2000000),include.lowest = T)

#
pos2<-rowMeans(cbind(aggregate(gend.sites[,2],by=list(pos.i,gend.sites[,1]),min,na.rm=T)[,3],
                     aggregate(gend.sites[,2],by=list(pos.i,gend.sites[,1]),max,na.rm=T)[,3])
)

atlasscan<-abs(finalhaps[,"X111"]-as.matrix(gend[,4]))/2
atlasscan<-aggregate(atlasscan,by=list(pos.i,gend.sites[,1]),mean,na.rm=T)
maisscan<-abs(finalhaps[,"X111"]-as.matrix(gend[,17]))/2
maisscan<-aggregate(maisscan,by=list(pos.i,gend.sites[,1]),mean,na.rm=T)

pdf("FIGURES_FILES_TABLES/lin1_ancestry.pdf",height=8,width=5)
par(mfrow=c(7,1), mai=c(.2,.5,0,.5))
for (i in unique(atlasscan[,2])){
  plot(1,ylim=c(-max(c(atlasscan[,3],maisscan[,3])),max(c(atlasscan[,3],maisscan[,3]))),xlim=c(0, max(pos2)),axes=F,xlab="",ylab="")
  points(pos2[atlasscan[,2]==i],atlasscan[atlasscan[,2]==i,3],type="h",col=viridis(20)[1])
  points(pos2[maisscan[,2]==i],-maisscan[maisscan[,2]==i,3],type="h",col=viridis(20)[10])
  axis(2,las=2)
  axis(4,las=2)
}
axis(1)
dev.off()

#Calculate two parents with the largets amount of near identical segments
outsimilar<-list()#matrix(NA,nrow=ncol(finalhaps),ncol=28)
outsimilar.v<-list()
for(j in 1:ncol(finalhaps)){
  print(j)
  regscan<-abs(finalhaps[,j]-as.matrix(gend[,1:28]))/2
  regscan2<-aggregate(regscan,by=list(pos.i,gend.sites[,1]),mean,na.rm=T)
  regscan3<-regscan2[,-(1:2)]<.01
  prim_p<-order(colSums(regscan3,na.rm=T),decreasing=T)[1]
  secon_p<-order(colSums(regscan3[!regscan3[,prim_p],],na.rm=T),decreasing=T)[1]
  outsimilar[[j]]<-colnames(regscan3)[c(prim_p,secon_p)]
  outsimilar.v[[j]]<-c(sum(regscan3[,prim_p],na.rm=T),sum(regscan3[!regscan3[,prim_p],secon_p],na.rm=T))
}
pands<-table(as.data.frame(do.call("rbind",outsimilar)))
pands_2<-melt(pands)
colnames(pands_2)<-c("Primary","Secondary","Count")
pands_2[,1]<-factor(pands_2[,1],levels=names(sort(rowSums(pands))))
pands_2[,2]<-factor(pands_2[,2],levels=names(sort(colSums(pands))))

#Calculate # of lineages with specific ancestry 
rowSums(pands)/sum(pands) # primary parent Fractions
colSums(pands)/sum(pands) # secondary parent Fractions
pands["Atlas","Arequipa"]/sum(pands) # 0.1954023
pands["Atlas","Maison_Carree_Carre_42"]/sum(pands) # 0.1494253


#Generate heatmap for most similar parents
pdf("FIGURES_FILES_TABLES/parent_ancestory.pdf")
heatmap <- ggplot(pands_2,
                  aes(Primary, Secondary, fill = Count)) +
  geom_tile(colour="white", size=1.5, stat="identity") + 
  scale_fill_distiller(palette = "Spectral",trans='sqrt') +
  geom_text(aes(label = round(Count, 2)), na.rm = TRUE) +
  labs(x = NULL, y = NULL) +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="white"),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(color="grey20", size=rel(1)),
    axis.text.y  = element_text(hjust=1),
    axis.text.x = element_text(angle=45,hjust=1,vjust = 1),
    legend.text = element_text(color="grey20", size=rel(1.3)),
    legend.background = element_rect(fill="white"),
    legend.position = "bottom",
    legend.title=element_blank()
  ) +
  coord_equal()
heatmap
p<-ggplot(data=pands_2, aes(x=Primary, y=Count)) +
  geom_bar(stat="identity", fill="grey") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="white"),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank()
  )
p
p<-ggplot(data=pands_2, aes(x=Secondary, y=Count)) +
  geom_bar(stat="identity", fill="grey") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="white"),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank()
  )
p
dev.off()
##################################################################################################################
##################################################################################################################
eigval<-read.table("PLINK/plink_string.eigenval")
eigvec<-read.table("PLINK/plink_string.eigenvec")
#
fam<-rep(NA,nrow(eigvec))
fam[1:28]<-0
fam[grep("^1_",eigvec[,1])]<-18
fam[grep("^2_",eigvec[,1])]<-28
fam[grep("^3_",eigvec[,1])]<-50
fam[grep("^7_",eigvec[,1])]<-58
#
pdf("FIGURES_FILES_TABLES/PCA_of_radseq.pdf",height=6,width=4)
par(mfrow=c(3,2),las=2)
plot(eigvec[,3:4],col="lightgrey",xlab="PC1",ylab="PC2",main="Parents",las=2)
points(eigvec[fam==0,3:4],pch=21,bg=brewer.pal(9,"Set1")[2])
points(eigvec[fam==0,3:4][4,],pch=21,bg=brewer.pal(9,"Set1")[1])

plot(eigvec[,3:4],col="lightgrey",xlab="PC1",ylab="PC2",main="F18",las=2)
points(eigvec[fam==18,3:4],pch=21,bg=brewer.pal(9,"Set1")[2])
plot(eigvec[,3:4],col="lightgrey",xlab="PC1",ylab="PC2",main="F28",las=2)
points(eigvec[fam==28,3:4],pch=21,bg=brewer.pal(9,"Set1")[2])
plot(eigvec[,3:4],col="lightgrey",xlab="PC1",ylab="PC2",main="F50",las=2)
points(eigvec[fam==50,3:4],pch=21,bg=brewer.pal(9,"Set1")[2])
plot(eigvec[,3:4],col="lightgrey",xlab="PC1",ylab="PC2",main="F58",las=2)
points(eigvec[fam==58,3:4],pch=21,bg=brewer.pal(9,"Set1")[2])
plot(eigvec[,3:4],col="lightgrey",xlim=c(-.013,-.01),ylim=c(-.001,.01),xlab="PC1",ylab="PC2",main="F58",las=2)
points(eigvec[fam==58,3:4],pch=21,bg=brewer.pal(9,"Set1")[2])
dev.off()
