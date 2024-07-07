setwd('~/bigdata/BARLEY_CCII_FULL/DATA/OUTPUT/')
###########################libraries################################################
library(ape)
####################################################################################
####################################################################################
#Build tree from distance matrix for parents
dm<-read.table("PLINK/plink_file_prune.dist")
colnames(dm)<-rownames(dm)<-read.table("PLINK/plink_file_prune.dist.id")[,1]
dmfull<-read.table("PLINK/plink_file_noprune.dist")

####################################################################################
colnames(dmfull)<-rownames(dmfull)<-read.table("PLINK/plink_file_noprune.dist.id")[,1]
dm.out<-dmfull/(13271163*2)
dm.out[upper.tri(dm.out,diag=F)]<-1-dm.out[upper.tri(dm.out,diag=F)]
write.table(round(dm.out,digits=3),"TABLE/dm_matrix.txt",quote=F)
#
pca<-read.table("PLINK/plink_file_prune.eigenvec")
varex<-read.table("PLINK/plink_file_prune.eigenval")[,1]
hybrids<-c("Flynn",
           "Meloy",
           "Alpha",
           "GoldenPheasant",
           "Glabron"
)

#Plot
pdf("FIGURE/Parental_relationships.pdf",width=10,height=10)
par(mfrow=c(2,2))
plot(nj(as.dist(dm[!colnames(dm)%in%hybrids,!rownames(dm)%in%hybrids])),type="u",lab4ut="axial")
barplot(varex/sum(varex),las=2,ylab="Fraction Variance Explained")
plot(pca[,3:4],xlab="PC1",ylab="PC2")
text(pca[,3:4],pca[,1])
plot(pca[,5:6],xlab="PC3",ylab="PC4")
text(pca[,5:6],pca[,1])
dev.off()
