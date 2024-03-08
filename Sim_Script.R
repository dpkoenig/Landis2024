library("data.table")
################################################################################
#function#
################################################################################
simrecomb<-function(x){
  recomb<-list()
  for (i in 1:length(x)){
    tparents<-sample(1:28,2,replace=F)
    breaks<-sort(sample(2:(nrow(x[[i]])-1),1,replace=F))
    recomb[[i]]<-c(x[[i]][1:breaks[1],tparents[1]],
                   x[[i]][(breaks[1]+1):nrow(x[[i]]),tparents[2]])  
  }
  return(unlist(recomb))
}
################################################################################
################################################################################
args = commandArgs(trailingOnly=TRUE)
print(length(args))
if (length(args)<6 | length(args)>6) {
  stop("Six arguments must be supplied (input file).n", call.=FALSE)
}
################################################################################
################################################################################
#1: Progeny AFS file
#2: Parent Genotypes file
#3: Founder Population Size file
#4: Site down sample size
#5: Rep #
#6: Out file
################################################################################
################################################################################
afs_o<-as.data.frame(fread(args[1]))
parents_o<-as.matrix(fread(args[2]))
fsize<-as.numeric(args[3])
ssize<-as.numeric(args[4])
repnum<-as.numeric(args[5])

#Generate output matrix
final_out<-matrix(NA,ncol=7,nrow=repnum)

#Begin simulations
for (j in 1:repnum){
  #site downsampling
  if (ssize>0){
    dsamp<-sample(1:nrow(afs_o),size=ssize,replace=F)
    dsamp<-sort(dsamp)
    afs<-afs_o[dsamp,]
    parents<-parents_o[dsamp,]
  }
  
  #Calculate observed data and split parent genotypes by chromosome
  poolc<-afs[,seq(6,12,2)]+afs[,seq(5,12,2)]
  posd<-split(afs[,2],afs[,1])
  gchro<-split.data.frame(parents,afs[,1])
  af<-afs[,seq(6,12,2)]/poolc
  he<-2*af*(1-af)
  
  #Simulate random evolution of lines and calculate the remaining diversity at
  #the apporpiate time points
  poprep<-list()
  poprep[[2]]<-1:fsize
  for (i in 3:58){
    poprep[[i]]<-sample(poprep[[i-1]],fsize,replace=T)
  }
  simlines<-unique(unlist(poprep[c(18,28,58)]))
  
  #For each line present in sampled generations, generate a random two parent recombinant
  #one recombination per chromosome
  outrecomb<-matrix(NA,ncol=length(simlines),nrow=nrow(parents))
  for (i in 1:ncol(outrecomb)){
    outrecomb[,i]<-simrecomb(gchro)
  }
  colnames(outrecomb)<-paste("X",unique(unlist(poprep[c(18,28,58)])),sep="")
  
  #Calculate He for each generation based on simulated individuals
  Heout<-c(NA,NA,NA)
  n<-1
  for (i in c(18,28,58)){
    lineabun<-table(poprep[[i]])
    names(lineabun)<-paste("X",names(lineabun),sep="")
    tmp<-outrecomb[,names(lineabun)]
    #
    tmp<-t(t(tmp)*as.vector(lineabun))
    tmpaf<-rowSums(tmp)/(2*sum(lineabun))
    Heout[n]<-mean(2*tmpaf*(1-tmpaf))
    n<-n+1
  }
  final_out[j,]<-c(colMeans(he),Heout)
}
write.table(final_out,args[6],col.names = F, row.names = F, quote=F)