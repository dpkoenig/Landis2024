#Accesory functions for afs analysis
###########################Function#################################################
#Calculate nei's genetic distance
neid<-function(x,y){
  Jx<-mean(x^2+(1-x)^2,na.rm=T)
  Jy<-mean(y^2+(1-y)^2,na.rm=T)
  Jxy<-mean(x*y+(1-x)*(1-y),na.rm=T)
  return(-log(Jxy/sqrt(Jx*Jy)))
}
#Set colors for plotting
getlinecols<-function(x){
  z<-read.table("../INPUT/PARENT_GROUPS.txt")
  #
  rownames(z)<-z[,1]
  #Convert to colors 
  z[z[,2]=="H" | z[,2]=="U",]<-brewer.pal(12,"Set3")[9]
  z[z[,2]=="E2" | z[,2]=="E6",]<-brewer.pal(12,"Set3")[5]
  z[z[,2]=="NoA",]<-brewer.pal(12,"Set3")[6]
  z[z[,2]=="A6",]<-brewer.pal(12,"Set3")[1]
  zz<-z[,2]
  names(zz)<-rownames(z)
  return(zz)
}
