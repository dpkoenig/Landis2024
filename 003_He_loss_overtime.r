#Script to explore and plot loss of genetic variation in CCII over time
#Older simulation approach shown below
##########################################################################
library(data.table)
#
setwd('~/bigdata/BARLEY_CCII_FULL/DATA/OUTPUT/')
##########################################################################
#Functions
##########################################################################
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
runone<-function(x,y,z){
  thisgen<-x
  outgen<-list()
  outgen[[2]]<-1:y
  for (i in 3:z){
    thisgen<-sample(thisgen,size=y,replace=T)
    outgen[[i]]<-thisgen
  }
  return(outgen)
}
sim_afs<-function(x){
  outsim<-table(factor(rep(1:28,round(p_size/28))[x],levels=1:28))
  simafs<-rowSums(gtnumsim*outsim)/sum(outsim*2)
  return(mean(2*simafs*(1-simafs)))
}
calcHE<-function(x){
  thecov<-x[,seq(1,ncol(x),2)]+x[,seq(2,ncol(x),2)]
  thefq<-x[,seq(2,ncol(x),2)]/thecov
  return(colMeans(2*thefq*(1-thefq)))
}
##########################################################################
#Data Read In
##########################################################################
ac<-fread("FILTERED/filtered_ac.txt", stringsAsFactors = F)
ac<-as.data.frame(ac,stingsAsFactors=F)
gtnum<-fread("FILTERED/filtered_gtnum.txt", stringsAsFactors = F)
gtnum<-as.data.frame(gtnum,stingsAsFactors=F)

##########################################################################
#Prep data and use only complete sites
##########################################################################
no.miss.i<-rowSums(is.na(gtnum))==0
acsim<-ac[no.miss.i,]
gtnumsim<-gtnum[no.miss.i,]
gtnumsim<-as.matrix(gtnumsim)/2
##########################################################################
#Observed He over time
##########################################################################
obsHe<-calcHE(acsim[,5:12])

##########################################################################
#linear model
##########################################################################
obsHemod<-data.frame(list(gen=c(2,18,28,58),He=obsHe))
Helm<-lm(He~gen,obsHemod)
summary(Helm)# p = 0.001220 slope = -0.0031822 R2 =  0.9963 
-coefficients(Helm)[1]/coefficients(Helm)[2]#70.58553

##########################################################################
#Run Simulations
##########################################################################
#Old Simulation approach no recombination
#sim10000<-list()
# for (i in 1:1000){
#   print(i)
#   p_size<-10000
#   startp<-1:p_size
#   thisrun<-runone(startp,p_size,58)
#   gtdown<-gtnumsim[sample(1:nrow(gtnumsim),size=100000,replace=F),]
#   outHe<-rep(NA,58)
#   for (j in c(18,28,58)){
#     outsim<-table(factor(rep(1:28,round(p_size/28))[thisrun[[j]]],levels=1:28))
#     allelefreq<-rowSums(sweep(gtdown, MARGIN=2, outsim, `*`))/p_size
#     outHe[j]<-mean(2*allelefreq*(1-allelefreq))
#   }
#   sim10000[[i]]<-outHe
# }
# 
# #
# 
# sim1000<-list()
# for (i in 1:1000){
#   print(i)
#   p_size<-1000
#   startp<-1:p_size
#   thisrun<-runone(startp,p_size,58)
#   gtdown<-gtnumsim[sample(1:nrow(gtnumsim),size=100000,replace=F),]
#   outHe<-rep(NA,58)
#   for (j in c(18,28,58)){
#     outsim<-table(factor(rep(1:28,round(p_size/28))[thisrun[[j]]],levels=1:28))
#     allelefreq<-rowSums(sweep(gtdown, MARGIN=2, outsim, `*`))/p_size
#     outHe[j]<-mean(2*allelefreq*(1-allelefreq))
#   }
#   sim1000[[i]]<-outHe
# }
# 
# #
# 
# sim100<-list()
# for (i in 1:1000){
#   print(i)
#   p_size<-100
#   startp<-1:p_size
#   thisrun<-runone(startp,p_size,58)
#   gtdown<-gtnumsim[sample(1:nrow(gtnumsim),size=100000,replace=F),]
#   outHe<-rep(NA,58)
#   for (j in c(18,28,58)){
#     outsim<-table(factor(rep(1:28,round(p_size/28))[thisrun[[j]]],levels=1:28))
#     allelefreq<-rowSums(sweep(gtdown, MARGIN=2, outsim, `*`))/p_size
#     outHe[j]<-mean(2*allelefreq*(1-allelefreq))
#   }
#   sim100[[i]]<-outHe
# }
# write.table(do.call("rbind",sim100),"TABLE/simulated_He_100.txt",quote=F,row.names = F,col.names = F)
# write.table(do.call("rbind",sim1000),"TABLE/simulated_He_1000.txt",quote=F,row.names = F,col.names = F)
# write.table(do.call("rbind",sim10000),"TABLE/simulated_He_10000.txt",quote=F,row.names = F,col.names = F)
# ##########################################################################
# #Plot Simulation Data
# ##########################################################################
# 
# pdf("FIGURE/Simulation_He.pdf")
# plot(0,
#      type="n",
#      ylim=c(0,.3),
#      xlim=c(0,80),
#      las=2,
#      xlab="Generation",
#      ylab="He",
#      axes=F,
#      xaxs = "i",
#      yaxs = "i")
# #polygon(c(0,0,58,58),c(0,.3,.3,0),col=brewer.pal(9,"Blues")[1],border = F)
# #polygon(c(58,58,80,80),c(0,.3,.3,0),col=brewer.pal(9,"Greys")[2],border = F)
# axis(1)
# axis(2,las=2)
# thecols<-viridis(3,.05)
# for(i in 1:length(sim100)){
#   lines(c(2,18,28,58),
#         c(obsHe[1],na.omit(sim100[[i]])),
#         col=thecols[1])
# }
# for(i in 1:length(sim1000)){
#   lines(c(2,18,28,58),
#         c(obsHe[1],na.omit(sim1000[[i]])),
#         col=thecols[2])
# }
# for(i in 1:length(sim10000)){
#   lines(c(2,18,28,58),
#         c(obsHe[1],na.omit(sim10000[[i]])),
#         col=thecols[3])
# }
# lines(c(2,18,28,58),obsHe,col="black",type="b",lwd=2)
# abline(lm(obsHe~c(2,18,28,58)),lty=2)
# dev.off()
#


#dsim<-do.call("rbind",sim100)
#sampbias<-function(x){mean(rbinom(x[1],1,prob=x[2]))}
#outsimd<-cbind(acsim.cov[,1],simafs)
#getbiased<-apply(outsimd,1,sampbias)
#mean(2*getbiased*(1-getbiased))
