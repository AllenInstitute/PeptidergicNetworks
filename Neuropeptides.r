# Figures - Smith et al


library(pheatmap)
library(RColorBrewer)
library(colorspace)
library(R.utils)
library(stringr)
library(superheat)
library(sm)
library(vioplot) 
 
# Read in processed data table for CPM and FPKM values from Tasic et al, Nature,  2018 Nov;563(7729):72-78 
# This table contains expression profiles for 23,823 cells from VISp and ALM cortical regions together with laminar location and cluster id.
# All figures and analyses can be generated from this code.

np_gpcr_cpm=as.data.frame(read.csv("/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/np_gpcr_cpm.csv"))
np_gpcr_cpm=as.data.frame(read.csv("F:\\Analysis\\Neuropeptides\\Submission\\np_gpcr_cpm.csv"))

# membership of cell types in Tasic and membership in major cell types
celltype=read.csv("/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/CellTypeAnnotation.csv",header=F)
celltype=read.csv("F:\\Analysis\\Neuropeptides\\Submission\\CellTypeAnnotation.csv",header=F)

 
# Gene classes of interest: NPP
nppm=c("Npy","Sst","Vip","Tac2", "Cck","Penk","Crh","Cort", "Tac1", "Pdyn", "Pthlh", "Pnoc","Trh","Grp","Rln1","Adcyap1","Nts","Nmb")
npph=toupper(nppm)

# NP-GPCR 
npgpcrm=c("Sstr2", "Npy2r", "Npy1r", "Grpr", "Cckbr", "Ntsr2", "Npy5r", "Nmbr", "Rxfp1", "Sstr4", 
          "Trhr","Sstr1","Adcyap1r1","Crhr1","Rxfp3","Oprl1","Crhr2","Tacr3","Oprk1","Tacr1",
          "Pth1r","Vipr1","Oprm1","Trhr2","Vipr2","Rxfp2","Oprd1","Ntsr1","Sstr3")


# twelve basic cell types for pooling results to coarser resolution
ctpool=c("IT","PT","NP","CT","L6b","Lamp5","Sncg","Serpinf1","VIP","Sst Chodl","Sst","Pvalb")
ctcolors=c("#A1D5BA","#7AB7AB","#BEDFCB","#9ECAE8","#91ACBA","#F2ABB7","#C387E2","#E6B8FF","#D5A8DF","#FEFF99","#E0B889","#EE98A4")
names(ctcolors)=ctpool

# totalnumber of types in taxonomy
cellind=1:133

# cell type membership data, types and colors for VISp and ALM regional analysis

allclust=read.csv("/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/All_Types.csv")
allclust=read.csv("F:\\Analysis\\Neuropeptides\\Submission\\All_Types.csv")
allclid=as.integer(allclust[,2])
alldend=as.integer(allclust[,1])

almclust=read.csv("/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/AlM_Types.csv")
almclust=read.csv("F:\\Analysis\\Neuropeptides\\Submission\\ALM_Types.csv")
almclid=as.integer(almclust[,2])
almdend=as.integer(almclust[,1])

visclust=read.csv("/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/VISp_Types.csv")
visclust=read.csv("F:\\Analysis\\Neuropeptides\\Submission\\VISp_Types.csv")
visclid=as.integer(visclust[,2])
visdend=as.integer(visclust[,1])
 
# Cell type lists for pooling to 12 major types, these give the index map for the 12 major cell types for each region

vislist=NULL
icount=0
for (i in 1:12) {
  itemp=intersect(which(celltype[,4]==i),visdend) 
  vislist[[i]]=c((icount+1):(icount+length(itemp)))
  icount=icount+length(itemp)
}

almlist=NULL
icount=0
for (i in 1:12) {
  itemp=intersect(which(celltype[,4]==i),almdend) 
  almlist[[i]]=c((icount+1):(icount+length(itemp)))
  icount=icount+length(itemp)
}


alllist=NULL
icount=0
for (i in 1:12) {
  itemp=intersect(which(celltype[,4]==i),alldend) 
  alllist[[i]]=c((icount+1):(icount+length(itemp)))
  icount=icount+length(itemp)
}

# table of 37 interaction pairs and 
inpairs=read.csv("F:\\Analysis\\Neuropeptides\\Submission\\CognatePairs.csv",header=FALSE)
inpairs=read.csv("/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/CognatePairs.csv",header=FALSE)



###################################################### 
#  Figure 1B. Fraction of Cells Expressing NPP 
###################################################### 

# Note: run all plotting subroutines below

# ALM - NPP
d1=subset(np_gpcr_cpm,brain_region_label=="ALM",nppm)
cvec=sum(apply(d1,1,sum)>0)
totcell=dim(d1)[1]
bimat=d1>0
almcnt=c(1,cvec/totcell,apply(bimat,2,sum)/totcell)
b1=almcnt[-1]
names(almcnt)[1]="ALM"
names(almcnt)[2]="Any NPP"
names(b1)[1]="Any"
pervec_alm=almcnt[c(-1,-2)]
 
 
# VISp - NPP
d1=subset(np_gpcr_cpm,brain_region_label=="VISp",nppm)
cvec=sum(apply(d1,1,sum)>0)
totcell=dim(d1)[1]
bimat=d1>0
almcnt=c(1,cvec/totcell,apply(bimat,2,sum)/totcell)
b2=almcnt[-1]
names(almcnt)[1]="VISp"
names(almcnt)[2]="Any NPP"
pervec_visp=almcnt[c(-1,-2)]

# dual plot
c1=rbind(b1,b2)
row.names(c1)=c("ALM","VISp")
cr1=terrain_hcl(19)[5]
cr2=terrain_hcl(19)[11]

par(las=2)
par(fig = c(0,1,0,1))
barplot(c1,col=rbind(rep(cr1,14),rep(cr2,14)),ylim=c(0,1),beside=TRUE,main="NPP Expressing Cells - ALM and VISp")
legend(40,0.8,legend=c("ALM","VISp"),col=c(cr1,cr2),pch=c(15,15),box.lty=0,pt.cex=1.5,cex=1.25)



########################################################## 
#  Figure 1C. Histogram for NPP Expressing Cells 
########################################################## 


d1=subset(np_gpcr_cpm,brain_region_label=="ALM",nppm)
totcell=dim(d1)[1]
dpos=d1>0
fvec=apply(dpos,1,sum)
tvec=table(fvec)
tvecn=tvec/sum(tvec[2:length(tvec)]) 
tvec1=tvecn


d1=subset(np_gpcr_cpm,brain_region_label=="VISp",nppm)
totcell=dim(d1)[1]
dpos=d1>0
fvec=apply(dpos,1,sum)
tvec=table(fvec)
tvecn=tvec/sum(tvec[2:length(tvec)]) 
tvec2=tvecn


# dual plot
c1=rbind(b1,b2)
row.names(c1)=c("ALM","VISp")
cr1=terrain_hcl(19)[5]
cr2=terrain_hcl(19)[11]

# in case one array is longer pad with terminal zeros
difl=length(tvec1)-length(tvec2)
if (difl > 0) {
  tvec2=c(tvec2,rep(0,difl))
  names(tvec2)=names(tvec1)
}
if (difl < 0) {
  tvec1=c(tvec1,rep(0,abs(difl)))
  names(tvec1)=names(tvec2)
}

outr=rbind(tvec1,tvec2)

par(las=1)
barplot(outr,col=rbind(rep(cr1,12),rep(cr2,12)),ylim=c(0,0.2),beside=TRUE,
        ylab="Histogram Frac. NPP expressing",cex.names=1.2)
title(main="Histogram of NPP expressing cells in ALM and VISp")
par(mar=c(5, 4, 4, 2.7) + 0.1)  
legend(30,0.15,legend=c("ALM","VISp"),col=c(cr1,cr2),pch=c(15,15),box.lty=0,pt.cex=1.5,cex=1.25) 


#########################################################
#  Figure 2B. Fraction of Cells Expressing NP-GPCR
######################################################### 


d1=subset(np_gpcr_cpm,brain_region_label=="ALM",npgpcrm)
cvec=sum(apply(d1,1,sum)>0)
totcell=dim(d1)[1]
bimat=d1>0
almcnt=c(1,cvec/totcell,apply(bimat,2,sum)/totcell)
b1=almcnt[-1]
names(almcnt)[1]="ALM"
names(b1)[1]="Any NP-GPCR"
pervec_alm=almcnt[c(-1,-2)]

d1=subset(np_gpcr_cpm,brain_region_label=="VISp",npgpcrm)
cvec=sum(apply(d1,1,sum)>0)
totcell=dim(d1)[1]
bimat=d1>0
almcnt=c(1,cvec/totcell,apply(bimat,2,sum)/totcell)
b2=almcnt[-1]
names(almcnt)[1]="ALM"
names(almcnt)[2]="Any NP-GPCR"
names(b1)[1]="Any NP-GPCR"
pervec_alm=almcnt[c(-1,-2)]


# dual plot
c1=rbind(b1,b2)
row.names(c1)=c("ALM","VISp")
cr1=rainbow_hcl(19)[5]
cr2=rainbow_hcl(19)[11]

par(las=2)
barplot(c1,col=rbind(rep(cr1,14),rep(cr2,14)),ylim=c(0,1),beside=TRUE,main="NP-GPCR Expressing Cells - ALM and VISp")
legend(70,0.8,legend=c("ALM","VISp"),col=c(cr1,cr2),pch=c(15,15),box.lty=0,pt.cex=1.5,cex=1.25)

###################################################################
#  Fraction of cell types expressing NP-GPCR and cognate pairs
###################################################################

commat=matrix(0,dim(inpairs)[1],118)

for (i in 1:dim(inpairs)[1]) {
  np1=as.character(inpairs[i,2])
  np2=as.character(inpairs[i,3])
  ipos=which(((np_gpcr_cpm[,np1]>0)*(np_gpcr_cpm[,np2]>0))==1)
  lipos=length(ipos)
  itab=table((np_gpcr_cpm[,"cluster_id"])[ipos])[alldend]
  commat[i,as.integer(names(itab))]=as.integer(itab)
} 

########################################################## 
#  Figure 2C. Histogram for NP-GPCR Expressing Cells 
########################################################## 


d1=subset(np_gpcr_cpm,brain_region_label=="ALM",npgpcrm)
totcell=dim(d1)[1]
dpos=d1>0
fvec=apply(dpos,1,sum)
tvec=table(fvec)
tvecn=tvec/sum(tvec[2:length(tvec)]) 
tvec1=tvecn
   
d1=subset(np_gpcr_cpm,brain_region_label=="VISp",npgpcrm)
totcell=dim(d1)[1]
dpos=d1>0
fvec=apply(dpos,1,sum)
tvec=table(fvec)
tvecn=tvec/sum(tvec[2:length(tvec)]) 
tvec2=tvecn
 
# dual plot
c1=rbind(tvec1,tvec2)
row.names(c1)=c("ALM","VISp")
cr1=rainbow_hcl(19)[5]
cr2=rainbow_hcl(19)[11]


# in case one array is longer pad with terminal zeros
difl=length(tvec1)-length(tvec2)
if (difl > 0) {
  tvec2=c(tvec2,rep(0,difl))
  names(tvec2)=names(tvec1)
}
if (difl < 0) {
  tvec1=c(tvec1,rep(0,abs(difl)))
  names(tvec1)=names(tvec2)
}

outr=rbind(tvec1,tvec2)



par(las=1)
barplot(outr,col=rbind(rep(cr1,12),rep(cr2,12)),ylim=c(0,0.2),beside=TRUE,
        ylab="Histogram Frac. NP-GPCR Expressing")
title(main="Histogram of NP-GPCR expressing cells in ALM and VISp")
par(mar=c(5, 4, 4, 2.7) + 0.1)  
legend(40,0.15,legend=c("ALM","VISp"),col=c(cr1,cr2),pch=c(15,15),box.lty=0,pt.cex=1.5,cex=1.25) 



#####################################################################################
#  Figure 3A. Basic expression plots and variability for NPP - Combined VISp and ALM
#####################################################################################


alltyp=133
clabel=paste("CT",c(1:alltyp),sep="")

nppmat=matrix (0,length(nppm),alltyp,dimnames=list(nppm,clabel))

for (i in 1:length(nppm)) {
  ip1=which(names(np_gpcr_cpm)==nppm[i])
  for (j in 1:alltyp) {
    ip2=which(np_gpcr_cpm[,6]==j)
     nppmat[i,j]=mean(np_gpcr_cpm[ip2,ip1],trim=0.05)
  }
} 


vmat=nppmat[,allclid]

# log before normalization
vmat=log10(vmat+1)


# before normalization and plotting, compute mean variation within major cell types compared to total variation

varmat=matrix(0,dim(vmat)[1],length(ctpool)+1,dimnames=list(rownames(vmat),c("All",ctpool)))

for (i in 1:dim(vmat)[1]) {
  varmat[i,1]=var(vmat[i,])/mean(vmat[i,])
}
for (i in 1:dim(vmat)[1]) {
  for (j in 1:length(ctpool)) {
    if (length(vmat[i,alllist[[j]]]) > 1) 
      if (mean(vmat[i,alllist[[j]]])>0) {
      varmat[i,j+1]=var(vmat[i,alllist[[j]]])/mean(vmat[i,alllist[[j]]])
      }
    }
} 

varmat_npp=varmat 
 

##################################################################### 
#  Figure 3C. Violin plots NPP variability over major types
#####################################################################

# create data frame for variability plots
varmatm=as.data.frame(varmat[,-11])
labn=colnames(varmatm)
longd=NULL
for (j in 1:dim(varmatm)[2]) {
  longd=rbind(longd,cbind(varmatm[,j],rep(labn[j],dim(varmatm)[1])))
}
longdf=as.data.frame(longd)
colnames(longdf)=c("dat","reg")
longdfp=longdf[as.numeric(as.character(longdf[,1]))>0,]


cvec=c(rgb(0.75,0.75,0.75),ctcolors)
names(cvec)[1]="All"
lcvec=NULL
for (i in 1:12) {
  lcvec=c(lcvec,rep(cvec[i],18))
}
longdf2=cbind(longdf,lcvec)


par(fig = c(0,1,0,1))
par(las=1)
vioplot(as.numeric(as.character(longdfp[longdfp$reg=="All",1])),as.numeric(as.character(longdfp[longdfp$reg=="IT",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="PT",1])),as.numeric(as.character(longdfp[longdfp$reg=="NP",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="CT",1])),as.numeric(as.character(longdfp[longdfp$reg=="L6b",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="Lamp5",1])),as.numeric(as.character(longdfp[longdfp$reg=="Sncg",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="Serpinf1",1])),as.numeric(as.character(longdfp[longdfp$reg=="VIP",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="Sst",1])),as.numeric(as.character(longdfp[longdfp$reg=="Pvalb",1])),col=cvec[-11],
        names=names(cvec)[-11],cex.axis=1.75,main="NPP expression variability",cex.main=1.5,axes=F,ylim=c(0,2.5),
        font.axis=1)
 points(2,max(varmatm[,2]),cex=1.25,pch=16)
 points(3,max(varmatm[,3]),cex=1.25,pch=16)
 points(4,max(varmatm[,4]),cex=1.25,pch=16)
 points(5,max(varmatm[,5]),cex=1.25,pch=16)
 points(6,max(varmatm[,6]),cex=1.25,pch=16)
 points(7,max(varmatm[,7]),cex=1.25,pch=16)
 points(8,max(varmatm[,8]),cex=1.25,pch=16)
 points(9,max(varmatm[,9]),cex=1.25,pch=16)
 points(10,max(varmatm[,10]),cex=1.25,pch=16)
 points(11,max(varmatm[,11]),cex=1.25,pch=16)
 points(11,max(varmatm[,11]),cex=1.25,pch=16)
 points(12,max(varmatm[,12]),cex=1.25,pch=16)
 points(12,max(varmatm[,12]),cex=1.25,pch=16)

 
 # barplot for ratios on same plot
 mbarv=vector("numeric",11)
 mall=mean(varmatm[,1])
 for (j in 2:dim(varmatm)[2]) {
   mbarv[j-1]=mean(varmatm[varmatm[,j]>0,j])/mall
 }
 
 par(fig = c(0.75,0.95,0.65,1.0),new=TRUE,las=2)
 barplot(mbarv,names.arg=colnames(varmatm)[-1],col=ctcolors[-10],cex.names=0.75,cex.axis=0.75,ylab="Relative CV")
 abline(h=mean(mbarv),lty=2)
 
 
# Supplementary data for Figure 3 variability plots
 
 
varmatn=varmat[,-c(1,11)]
sumvec=apply(varmatn,2,mean)
ov=order(sumvec,decreasing=T)
pmat=varmatn[,ov]
 
covo=names(sumvec[ov])
pmatn=pmat
mypal=colorRampPalette(c(rgb(1,1,1),rgb(0.25,0.25,0.25)),bias=3)(10)
mypal=rev(colorRampPalette(brewer.pal(8,"Spectral"),bias=1)(8))

c1=brewer.pal(8,"Spectral")[1]
c8=brewer.pal(8,"Spectral")[8]
mypal=rev(colorRampPalette(c(rgb(0.75,0.4,1,0.75), rgb(0.75,0.4,1,0)), alpha = TRUE)(8))

 
imname="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/cell_type_variability_npp.tiff"
imname="F:\\Analysis\\Neuropeptides\\Submission\\cell_type_variability_npp.tiff"

tiff(imname,width=1500,height=1100)
par(cex.axis=1.5)
superheat(X=pmatn,left.label.col=c("ivory","ivory2"),bottom.label.col=ctcolors[covo],legend.text.size=14,
          left.label.text.alignment="center",X.text=round(pmatn,2),X.text.size=8,heat.pal=c("white","tan1"),yt=sumvec[covo],
          yt.axis.name="Mean cell type CV / Global CV",yt.plot.type="bar",yt.axis.size=18,yt.axis.name.size=16,yt.obs.col=ctcolors[covo],
          left.label.text.size=8,bottom.label.text.size=9)
dev.off()


# normalalize to row maximum
for (i in 1:dim(vmat)[1])  {
  print(paste(rownames(nppmat)[i],round(max(vmat[i,]),3)))
  vmat[i,]=vmat[i,]/max(vmat[i,])
}


# run plotting subroutines first
plot_All(vmat,"NPP expression in ALM and VISp",1) 

# expression pooled to 12 major types
plot_All_pooled(vmat,"NPP expression in ALM and VISp",1)


############################################################
#  Figure 3B. Basic expression plots for NP-GPCR - Combined
############################################################ 

alltyp=133
clabel=paste("CT",c(1:alltyp),sep="")

nppmat=matrix(0,length(npgpcrm),alltyp,dimnames=list(npgpcrm,clabel))

for (i in 1:length(npgpcrm)) {
  ip1=which(names(np_gpcr_cpm)==npgpcrm[i])
  for (j in 1:alltyp) {
    ip2=which(np_gpcr_cpm[,6]==j)
    nppmat[i,j]=mean(np_gpcr_cpm[ip2,ip1],trim=0.05)
  }
} 


vmat=nppmat[,allclid]


# log before normalization
vmat=log10(vmat+1)


# before normalization and plotting, compute mean variation within major cell types compared to total variation

varmat=matrix(0,dim(vmat)[1],length(ctpool)+1,dimnames=list(rownames(vmat),c("All",ctpool)))

for (i in 1:dim(vmat)[1]) {
  varmat[i,1]=var(vmat[i,])/mean(vmat[i,])
}
for (i in 1:dim(vmat)[1]) {
  for (j in 1:length(ctpool)) {
    if (length(vmat[i,alllist[[j]]]) > 1) 
      if (mean(vmat[i,alllist[[j]]])>0) {
        varmat[i,j+1]=var(vmat[i,alllist[[j]]])/mean(vmat[i,alllist[[j]]])
      }
  }
}


varmat_npgpcr=varmat 
 

##################################################################### 
#  Figure 3D. Violin plots for NP-GPCR variability over major types
#####################################################################

varmatm=as.data.frame(varmat[,-11])
labn=colnames(varmatm)
longd=NULL
for (j in 1:dim(varmatm)[2]) {
  longd=rbind(longd,cbind(varmatm[,j],rep(labn[j],dim(varmatm)[1])))
}
longdf=as.data.frame(longd)

colnames(longdf)=c("dat","reg")
longdfp=longdf[as.numeric(as.character(longdf[,1]))>0,]


# violin plot for cv variability

cvec=c(rgb(0.75,0.75,0.75),ctcolors)
names(cvec)[1]="All"
lcvec=NULL
for (i in 1:12) {
  lcvec=c(lcvec,rep(cvec[i],29))
}
longdf2=cbind(longdf,lcvec)
 

par(fig = c(0,1,0,1))
par(las=1)
vioplot(as.numeric(as.character(longdfp[longdfp$reg=="All",1])),as.numeric(as.character(longdfp[longdfp$reg=="IT",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="PT",1])),as.numeric(as.character(longdfp[longdfp$reg=="NP",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="CT",1])),as.numeric(as.character(longdfp[longdfp$reg=="L6b",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="Lamp5",1])),as.numeric(as.character(longdfp[longdfp$reg=="Sncg",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="Serpinf1",1])),as.numeric(as.character(longdfp[longdfp$reg=="VIP",1])),
        as.numeric(as.character(longdfp[longdfp$reg=="Sst",1])),as.numeric(as.character(longdfp[longdfp$reg=="Pvalb",1])),col=cvec[-11],
        names=names(cvec)[-11],cex.axis=1.75,main="NP-GPCR expression variability",cex.main=1.5,axes=F,ylim=c(0,2),
        font.axis=1)
points(2,1.7113,cex=1.25,pch=16)
points(3,1.027,cex=1.55,pch=16)
points(4,0.821417,cex=1.55,pch=16)
points(5,0.7069,cex=1.25,pch=16)
points(6,0.95626,cex=1.25,pch=16)
points(7,1.4198,cex=1.25,pch=16)
points(7,1.1718,cex=1.25,pch=16)
points(8,0.9220,cex=1.25,pch=16)
points(9,1.2237,cex=1.25,pch=16)
points(10,1.127,cex=1.25,pch=16)
points(11,1.38148,cex=1.25,pch=16)
points(11,1.228,cex=1.25,pch=16)
points(12,1.446,cex=1.25,pch=16)
points(12,1.338,cex=1.25,pch=16)


mbarv=vector("numeric",11)
mall=mean(varmatm[,1])
for (j in 2:dim(varmatm)[2]) {
  mbarv[j-1]=mean(varmatm[varmatm[,j]>0,j])/mall
}

par(fig = c(0.75,0.95,0.65,1.0),new=TRUE,las=2)
barplot(mbarv,names.arg=colnames(varmatm)[-1],col=ctcolors[-10],cex.names=0.75,cex.axis=0.75,ylab="Relative CV")
abline(h=mean(mbarv),lty=2)


# Supplementary data for Figure 3 variability plots


varmatn=varmat[,-c(1,11)]
sumvec=apply(varmatn,2,mean)
ov=order(sumvec,decreasing=T)
pmat=varmatn[,ov]

covo=names(sumvec[ov])
pmatn=pmat
mypal=colorRampPalette(c(rgb(1,1,1),rgb(0.25,0.25,0.25)),bias=3)(10)
mypal=rev(colorRampPalette(brewer.pal(8,"Spectral"),bias=1)(8))

c1=brewer.pal(8,"Spectral")[1]
c8=brewer.pal(8,"Spectral")[8]
mypal=rev(colorRampPalette(c(rgb(0.75,0.4,1,0.75), rgb(0.75,0.4,1,0)), alpha = TRUE)(8))
 

imname="/volumes/PHILIPS UFD/Analysis/Neuropeptides/cell_type_variability_np-gpcr.tiff"
imname="F:\\Analysis\\Neuropeptides\\Submission\\cell_type_variability_np-gpcr.tiff"
 
tiff(imname,width=1500,height=1500)
par(cex.axis=1.5)
superheat(X=pmatn,left.label.col=c("ivory","ivory2"),bottom.label.col=ctcolors[covo],legend.text.size=14,
          left.label.text.alignment="center",X.text=round(pmatn,2),X.text.size=8,heat.pal=c("white","tan1"),yt=sumvec[covo],
          yt.axis.name="Mean cell type CV / Global CV",yt.plot.type="bar",yt.axis.size=18,yt.axis.name.size=16,yt.obs.col=ctcolors[covo],
          left.label.text.size=8,bottom.label.text.size=9)
dev.off()


for (i in 1:dim(vmat)[1])  {
  print(paste(rownames(nppmat)[i],round(max(vmat[i,]),3)))
  vmat[i,]=vmat[i,]/max(vmat[i,])
}

plot_All(vmat,"NP-GPCR expression in ALM and VISp",1) 
plot_All_pooled(vmat,"NP-GPCR expression in ALM and VISp",1)


colnames(vmat)=paste("CT",alldend,sep="")
write.csv(vmat,"F:\\Analysis\\Neuropeptides\\FigureData\\Figure3B.csv")


#  try plotting one large

combmat=rbind(varmat_npp,varmat_npgpcr)
combn=combmat[,-c(1,11)]

pheatmap(combn)

##################################
# 5A.  NPP: VISp Expression
################################## 
 

d1=subset(np_gpcr_cpm,brain_region_label=="VISp")

alltyp=133
clabel=paste("CT",c(1:alltyp),sep="")

nppmat=matrix(0,length(nppm),alltyp,dimnames=list(nppm,clabel))

for (i in 1:length(nppm)) {
  ip1=which(names(d1)==nppm[i])
  for (j in 1:alltyp) {
    ip2=which(d1[,6]==j)
    nppmat[i,j]=mean(d1[ip2,ip1],trim=0.05)
  }
} 

vmat=nppmat[,visclid]


vmat_npp_visp=vmat

# log before normalization
vmat=log10(vmat+1)


for (i in 1:dim(vmat)[1])  {
  print(paste(rownames(nppmat)[i],round(max(vmat[i,]),3)))
  vmat[i,]=vmat[i,]/max(vmat[i,])
}

plot_VISp(vmat,"NPP expression in VISp",1) 

colnames(vmat)=paste("CT",visdend,sep="")
write.csv(vmat,"F:\\Analysis\\Neuropeptides\\FigureData\\Figure5A.csv")


######################################
# 5B.  NP-GPCR: VISp Expression
###################################### 

d1=subset(np_gpcr_cpm,brain_region_label=="VISp")

alltyp=133
clabel=paste("CT",c(1:alltyp),sep="")

nppmat=matrix(0,length(npgpcrm),alltyp,dimnames=list(npgpcrm,clabel))


for (i in 1:length(npgpcrm)) {
  ip1=which(names(d1)==npgpcrm[i])
  for (j in 1:alltyp) {
    ip2=which(d1[,6]==j)
    nppmat[i,j]=mean(d1[ip2,ip1],trim=0.05)
  }
} 

vmat=nppmat[,visclid]
vmat_npgpcr_visp=vmat



vmat_npgpcr_visp=vmat

# log before normalization
vmat=log10(vmat+1)


for (i in 1:dim(vmat)[1])  {
  print(paste(rownames(nppmat)[i],round(max(vmat[i,]),3)))
  vmat[i,]=vmat[i,]/max(vmat[i,])
}

plot_VISp(vmat,"NP-GPCR expression in VISp",1) 


colnames(vmat)=paste("CT",visdend,sep="")
write.csv(vmat,"F:\\Analysis\\Neuropeptides\\FigureData\\Figure5A.csv")




##################################
# 5C.  NPP: ALM Expression
##################################


d1=subset(np_gpcr_cpm,brain_region_label=="ALM")

alltyp=133
clabel=paste("CT",c(1:alltyp),sep="")

nppmat=matrix(0,length(nppm),alltyp,dimnames=list(nppm,clabel))

for (i in 1:length(nppm)) {
  ip1=which(names(d1)==nppm[i])
  for (j in 1:alltyp) {
    ip2=which(d1[,6]==j)
    nppmat[i,j]=mean(d1[ip2,ip1],trim=0.05)
  }
} 

vmat=nppmat[,almclid]


vmat_npp_alm=vmat

# log before normalization
vmat=log10(vmat+1)

for (i in 1:dim(vmat)[1])  {
  print(paste(rownames(nppmat)[i],round(max(vmat[i,]),3)))
  vmat[i,]=vmat[i,]/max(vmat[i,])
}

plot_ALM(vmat,"NPP expression in ALM",1) 


colnames(vmat)=paste("CT",almdend,sep="")
write.csv(vmat,"F:\\Analysis\\Neuropeptides\\FigureData\\Figure5C.csv")


######################################
# 5D.  NP-GPCR: ALM Expression
######################################



alltyp=133
clabel=paste("CT",c(1:alltyp),sep="")

nppmat=matrix(0,length(npgpcrm),alltyp,dimnames=list(npgpcrm,clabel))


for (i in 1:length(npgpcrm)) {
  ip1=which(names(d1)==npgpcrm[i])
  for (j in 1:alltyp) {
    ip2=which(d1[,6]==j)
    nppmat[i,j]=mean(d1[ip2,ip1],trim=0.05)
  }
} 

vmat=nppmat[,almclid]
vmat_npgpcr_alm=vmat


vmat_npgpcr_alm=vmat

# log before normalization
vmat=log10(vmat+1)


for (i in 1:dim(vmat)[1])  {
  print(paste(rownames(nppmat)[i],round(max(vmat[i,]),3)))
  vmat[i,]=vmat[i,]/max(vmat[i,])
}


plot_ALM(vmat,"NP-GPCR expression in ALM",1) 


colnames(vmat)=paste("CT",almdend,sep="")
write.csv(vmat,"F:\\Analysis\\Neuropeptides\\FigureData\\Figure5D.csv")


# Calculate intersection of cell types correlation - NPP
commont=intersect(colnames(vmat_npp_alm),colnames(vmat_npp_visp))
amat=as.vector(vmat_npp_alm[,commont])
bmat=as.vector(vmat_npp_visp[,commont])
plot(amat,bmat,pch=16,cex=0.5)
cor(amat,bmat)


# Calculate intersection of cell types correlation - NP-GPCR
commont=intersect(colnames(vmat_npgpcr_alm),colnames(vmat_npgpcr_visp))
amat=as.vector(vmat_npgpcr_alm[,commont])
bmat=as.vector(vmat_npgpcr_visp[,commont])

par(fig = c(0,1,0,1))
plot(log(amat+1),log(bmat+1),pch=16,cex=0.5,axes=F,xlab="ALM",ylab="VISp",main="Expression in common ALM and VISp cell types")
axis(1)
axis(2)
cor.test(amat,bmat)
amod=lm(log(bmat+1)~log(amat+1))
abline(amod$coefficients,col=2,lw=2)
rsq=summary(amod)$adj.r.squared
pval="p<2e-16"
text(2,5.75,paste("Adjusted R-squared = ",round(rsq,3)))
text(1.5,5.5,"p<2e-16")
 
    

#############################################################################################################
# # Figure 6. Expression profiles of conjugate NPP/NP-GPCR pairs predict neuropeptidergic signaling networks
#############################################################################################################
 
maxvec_alm=vector("numeric",dim(inpairs)[1])
maxvec_visp=vector("numeric",dim(inpairs)[1])

# assumes vat_npp_visp etc. matrices above

regionvec=c("VISp","ALM") 

for (region in regionvec) {

 for (ival in 1:dim(inpairs)[1]) {
  
  print (paste("Printing",region,ival))
  
  jval = as.character(inpairs[ival,2])
  kval = as.character(inpairs[ival,3])
  
  if (region == "VISp") {
     outerprod=outer(vmat_npp_visp[jval,],vmat_npgpcr_visp[kval,])
     outerprod=log10(outerprod+1)
     maxvec_visp[ival]=max(outerprod)
  }
  
  if (region == "ALM")  {
     outerprod=outer(vmat_npp_alm[jval,],vmat_npgpcr_alm[kval,])
     outerprod=log10(outerprod+1)
     maxvec_alm[ival]=max(outerprod)
  }
  
  # make title
  normv=round(max(outerprod),2)
  icode=paste("(",region,", ",as.character(inpairs[ival,4]),", ",normv,")",sep="")
  inum=paste(ival,":",sep="")
  
  ititle=paste(inum,jval,"->",kval, "  ", icode,sep=" ")
  
  fdir="F:\\Analysis\\Neuropeptides\\Submission\\Interaction Plots\\"
 # fdir="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Interaction Plots/"
  
  fname=paste(fdir,region,"_",ival,"_",as.character(inpairs[ival,1]),".tiff",sep="")
  
  tiff(fname,width=1000,height=1000)
  
  if (region == "VISp") {plot_Interaction_VISp(outerprod,ititle,gamma=1)}
  if (region == "ALM")  {plot_Interaction_ALM(outerprod,ititle,gamma=1)}
  
  dev.off()
  
  f2dir="F:\\Analysis\\Neuropeptides\\Submission\\Interaction Matrices\\"
  #f2dir="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Interaction Matrices/"
 
   f2name=paste(f2dir,ival,"_",region,"_",inpairs[ival,1],".csv",sep="")
  write.table(outerprod,file=f2name,col.names=F,row.names=F,quote=F,sep=",")
      
  
 }
}
  
 
# write out maximum scale normalization 
df=data.frame(cbind(as.character(inpairs[,1]),round(maxvec_alm,3),round(maxvec_visp,3)))
write.csv(df,"F:\\Analysis\\Neuropeptides\\Submission\\CognatePair_interaction_normalization.csv")
#write.csv(df,"/volumes/PHILIPS UFD/Analysis/Neuropeptides/Normalization.csv")



############################################
# Figure 6e: Pooled cell type representation
############################################

indir = "F:\\Analysis\\Neuropeptides\\Submission\\Interaction Matrices\\"
infiles=list.files(indir)

poolmat=matrix(0,12,12,dimnames=list(ctpool,ctpool))

plist_ALM=NULL
plist_VISp=NULL
ialm=0
ivisp=0

for (i in 1:length(infiles)) {
  
  print (paste("Printing Pooled",region,i))
  
  idf=read.csv(paste(indir,infiles[i],sep=""),header=F)
  imat=as.matrix(idf)
  
  fn1=unlist(strsplit(infiles[i],"[.]"))[1]
  ft=unlist(strsplit(fn1,"_"))
  fnum=as.integer(ft[1])
  region=as.character(ft[2])
  ipair=unlist(strsplit(as.character(ft[3]),"-"))  # construct interaction pair and read off type
  iloc=which(as.character(inpairs[,1])==as.character(ft[3]))
  itype=as.character(inpairs[iloc,4])
  
  for (j in 1:length(ctpool)) {
    for (k in 1:length(ctpool)) {
      
      if (region == "ALM")  {suml=almlist}
      if (region == "VISp") {suml=vislist}
    
      poolmat[j,k] = mean(imat[suml[[j]],suml[[k]]])
    }
  }
  
  if (region == "ALM")  {
    ialm=ialm+1
    plist_ALM[[ialm]]=poolmat
  }
  if (region == "VISp") {
    ivisp=ivisp+1
    plist_VISp[[ivisp]]=poolmat
  }
  

  maxv=round(max(poolmat),2)  # normalization factor
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(ctpool)
  )
  rownames(annotation_col)= ctpool
  
  ann_colors = list(CellType = ctcolors)
   
   
  mpal = colorRampPalette(brewer.pal(8,"PuBuGn"),bias=bval)(6)
   
  # construct title 
  titp = paste(fnum,": ",ipair[1],"->",ipair[2], "  (",region,", ",itype,", ",maxv, ")",sep="")
  fdir ="F:\\Analysis\\Neuropeptides\\Submission\\Interaction Plots Pooled\\"
  fname=paste(fdir,region,"_",fnum,"_",ipair[1],"-",ipair[2],"_Pooled.tiff",sep="")
  
  tiff(fname,width=800,height=800)
   
  
  pheatmap(poolmat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=20,fontsize_col=20,fontsize=20,border_color=9,color=rev(mpal),
           main=titp,annotation_col=annotation_col,annotation_row=annotation_col,annotation_colors = ann_colors,annotation_legend=F,
           show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=F,annotation_names_col=F,legend=F)
  
  dev.off()
  
  
}
 


##############################################################
#  Figure 6 (Supplmental) Clustering of interactions matrices
##############################################################

# make distance matrix between interaction matrices


fpath="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Interaction Matrices/"
fpath="F:\\Analysis\\Neuropeptides\\Interaction Matrices\\"

files=list.files(fpath)

matlistALM=NULL
matlistVISp=NULL
macount=0
mvcount=0
nameALM=vector("character",37)
nameVISp=vector("character",37)


corALM=matrix(0,37,37,dimnames=list(inpairs[,1],inpairs[,1]))
corVISp=matrix(0,37,37,dimnames=list(inpairs[,1],inpairs[,1]))

for (i in 1:length(files)) {
  
  ireg=unlist(strsplit(unlist(strsplit(files[i],"-"))[1],"_"))[2]
  inpp=unlist(strsplit(unlist(strsplit(files[i],"-"))[1],"_"))[3]
  igpcr=unlist(strsplit(unlist(strsplit(files[i],"-"))[2],"[.]"))[1] 
  ipos=which ((as.character(inpairs[,2])==inpp) & (as.character(inpairs[,3])==igpcr))
  icode=as.character(inpairs[ipos,4])
  ival=paste(inpp,igpcr,sep="-")
  
  if (ireg=="ALM")  {
    macount=macount+1
    matlistALM[[macount]]=as.matrix(read.csv(paste(fpath,files[i],sep=""),header=F))
    nameALM[macount]=ival
  }
  if (ireg=="VISp")  {
    mvcount=mvcount+1
    matlistVISp[[mvcount]]=as.matrix(read.csv(paste(fpath,files[i],sep=""),header=F))
    nameVISp[mvcount]=ival
  }
}

for (i in 1:macount)  {
  for (j in i:macount) {
    corALM[nameALM[i],nameALM[j]]=cor(as.vector(matlistALM[[i]]),as.vector(matlistALM[[j]]))
  }
}

for (i in 1:macount)  {
  for (j in i:macount) {
    corALM[j,i]=corALM[i,j]
    corALM[i,i]=0.0
  }
}



imname="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/ALM_interactions.tiff"
tiff(imname,width=1500,height=1500) 
superheat(X=corALM,row.dendrogram=T,col.dendrogram=T,bottom.label.text.angle=90,left.label.text.size=6,
          legend.text.size=16,bottom.label.text.size=6,title="ALM interactions",title.size=14)
dev.off()


for (i in 1:mvcount)  {
  for (j in i:mvcount) {
    corVISp[nameVISp[i],nameVISp[j]]=cor(as.vector(matlistVISp[[i]]),as.vector(matlistVISp[[j]]))
  }
}

for (i in 1:mvcount)  {
  for (j in i:mvcount) {
    corVISp[j,i]=corVISp[i,j]
    corVISp[i,i]=0.0
  }
}



imname="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/VISp_interactions.tiff"
tiff(imname,width=1500,height=1500) 
superheat(X=corVISp,row.dendrogram=T,col.dendrogram=T,bottom.label.text.angle=90,left.label.text.size=6,
          legend.text.size=16,bottom.label.text.size=6,title="VISp interactions",title.size=14)
dev.off()


s=svd(corALM)
eigval=s$d
alm_var=cumsum(eigval*eigval)/sum(eigval*eigval)
barplot(eigval*eigval,main="ALM ")


s=svd(corVISp)
eigval=s$d
visp_var=cumsum(eigval*eigval)/sum(eigval*eigval)



###################################################################################
# Figure 6 (Suppplemental): localization of interaction matrices to major cell types
###################################################################################


fpath="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Interaction Matrices/"
fpath="F:\\Analysis\\Neuropeptides\\Submission\\Interaction Matrices\\"


files=list.files(fpath)


locMat=matrix(0,37,2)
colnames(locMat)=c("ALM","VISp")
rownames(locMat)=inpairs[,1]

for (i in 1:length(files)) {
  
  ireg=unlist(strsplit(unlist(strsplit(files[i],"-"))[1],"_"))[2]
  inpp=unlist(strsplit(unlist(strsplit(files[i],"-"))[1],"_"))[3]
  igpcr=unlist(strsplit(unlist(strsplit(files[i],"-"))[2],"[.]"))[1] 
  ipos=which ((as.character(inpairs[,2])==inpp) & (as.character(inpairs[,3])==igpcr))
  icode=as.character(inpairs[ipos,4])
  ival=paste(inpp,igpcr,sep="-")
  
  
  fin=as.matrix(read.csv(paste(fpath,files[i],sep=""),header=F))
  
  vdat=as.vector(fin)
  hist(vdat,nclass=50,col=4)
  
  qthresh=quantile(vdat[vdat>0],probs=seq(0,1,0.1))[10]  
  
  # pool data exceeding threshhold to major cell types 
  
  #plot_Interaction_ALM(fin,"pre-thresh",1)
  
  for (j in 1:dim(fin)[1])  {
    for (k in 1:dim(fin)[2])  {
      if (fin[j,k]<qthresh)  {fin[j,k]=0.0}
    }
  }
  
  #plot_Interaction_ALM(fin,"post-thresh",1)
  
  omat=poolMat(fin,ireg)
  
  pfill=length(as.vector(omat[omat>0.0]))/144
  pfill2=length(as.vector(fin[fin>0.0]))*100.0/(dim(fin)[1]*dim(fin)[2])
  
  print(c(ireg,ival,pfill,pfill2))
  
  if (ireg=="ALM")  {locMat[ival,1]=round(pfill,3)*100.0}
  if (ireg=="VISp")  {locMat[ival,2]=round(pfill,3)*100.0}
  
}


# compute correlation of pooled interaction matrices between ALM and VISp
intercor=vector("numeric",dim(inpairs)[1])
for (i in 1:dim(inpairs)[1]) {
  intercor[i]=cor(as.vector(plist_ALM[[i]]),as.vector(plist_VISp[[i]]))
}

odat=cbind(locMat,intercor)
write.csv(odat,"F:\\Analysis\\Neuropeptides\\Submission\\Localization.csv")

mvals=apply(odats,2,max)
nodats=odats
nodats[,1]=mvals[1]-odats[,1] 
nodats[,2]=mvals[2]-odats[,2] 

par(las=2)
par(mar=c(7, 5, 4, 3.5) + 0.1)
odats=odat[order(odat[,3],decreasing = T),]
barplot(t(nodats[,1:2]),beside=T,horiz=F,ylab="Localization",xlab="",col=c("#fdbb84","#43a2ca"),cex.names=0.8,ylim=c(0,40),
        cex.main=2,main="Localization scores and Correlation: ALM and VISp")
legend(locator(1),legend=c("ALM","VISp"),col=c("#fdbb84","#43a2ca"),pch=15,bty="n",cex=1.5)
par(new=T)
plot(nodats[,3],axes=F,type="l",lwd=4,col=3,xlab="",ylab="")
axis(4)





#############################################################
# Interaction Classes with I,S,Q variables
#############################################################

# read in files by region maintaining counts exceeeding threshold by galpha (Q, I, S) type

fpath="F:\\Analysis\\Neuropeptides\\Interaction Matrices\\"
fpath="/volumes/PHILIPS UFD/Analysis/Neuropeptides/Submission/Interaction Matrices/"

files=list.files(fpath)


Smat_ALM=matrix(0,length(almclid),length(almclid))
Qmat_ALM=matrix(0,length(almclid),length(almclid))
Imat_ALM=matrix(0,length(almclid),length(almclid))


Smat_VISp=matrix(0,length(visclid),length(visclid))
Qmat_VISp=matrix(0,length(visclid),length(visclid))
Imat_VISp=matrix(0,length(visclid),length(visclid))


icnt=0; scnt=0; qcnt=0;

for (i in 1:length(files))  {
  
  print(paste("Computing I,S,Q",i))
  
  ireg=unlist(strsplit(unlist(strsplit(files[i],"-"))[1],"_"))[2]
  inpp=unlist(strsplit(unlist(strsplit(files[i],"-"))[1],"_"))[3]
  igpcr=unlist(strsplit(unlist(strsplit(files[i],"-"))[2],"[.]"))[1] 
  ipos=which ((as.character(inpairs[,2])==inpp) & (as.character(inpairs[,3])==igpcr))
  icode=as.character(inpairs[ipos,4])

  fin=read.csv(paste(fpath,files[i],sep=""),header=F)
  
  if (ireg ==  "ALM") {
    if (icode=="S") {
      Smat_ALM=Smat_ALM+10^fin
      scnt=scnt+1
     }
    if (icode=="Q") {
      Qmat_ALM=Qmat_ALM+10^fin
      qcnt=qcnt+1
      }
    if (icode=="I") {
      Imat_ALM=Imat_ALM+10^fin
      icnt=icnt+1
    }
  }
  
  if (ireg ==  "VISp") {
    if (icode=="S") {
      Smat_VISp=Smat_VISp+10^fin
      scnt=scnt+1
    }
    if (icode=="Q") {
      Qmat_VISp=Qmat_VISp+10^fin
      qcnt=qcnt+1
    }
    if (icode=="I") {
      Imat_VISp=Imat_VISp+10^fin
      icnt=icnt+1
    }
  }
}

icnt=icnt/2
scnt=scnt/2
qcnt=qcnt/2

# mean signal by I,S,Q

 Imat_ALM=Imat_ALM/icnt  
 Smat_ALM=Smat_ALM/scnt
 Qmat_ALM=Qmat_ALM/qcnt
  
 Imat_VISp=Imat_VISp/icnt
 Smat_VISp=Smat_VISp/scnt
 Qmat_VISp=Qmat_VISp/qcnt
 
 
 regionvec=c("VISp","ALM") 
 
 
 if (region == "ALM") {
   nms=paste("CT",almdend, sep="")
   rownames(Imat_ALM)=nms
   colnames(Imat_ALM)=nms
   rownames(Smat_ALM)=nms
   colnames(Smat_ALM)=nms
   rownames(Qmat_ALM)=nms
   colnames(Qmat_ALM)=nms
   
   write.csv(log10(Imat_ALM+1),paste("F:\\Analysis\\Neuropeptides\\ISQ\\Imat_ALM.csv",sep=""))  
   write.csv(log10(Smat_ALM+1),paste("F:\\Analysis\\Neuropeptides\\ISQ\\Smat_ALM.csv",sep=""))
   write.csv(log10(Qmat_ALM+1),paste("F:\\Analysis\\Neuropeptides\\ISQ\\Qmat_ALM.csv",sep=""))
 }
 
 
 
 if (region == "VISp") {
   nms=paste("CT",visdend, sep="")
   rownames(Imat_VISp)=nms
   colnames(Imat_VISp)=nms
   rownames(Smat_VISp)=nms
   colnames(Smat_VISp)=nms
   rownames(Qmat_VISp)=nms
   colnames(Qmat_VISp)=nms

   write.csv(log10(Imat_VISp+1),paste("F:\\Analysis\\Neuropeptides\\ISQ\\Imat_VISp.csv",sep=""))
   write.csv(log10(Smat_VISp+1),paste("F:\\Analysis\\Neuropeptides\\ISQ\\Smat_VISp.csv",sep=""))
   write.csv(log10(Qmat_VISp+1),paste("F:\\Analysis\\Neuropeptides\\ISQ\\Qmat_VISp.csv",sep=""))
}

# Full Size ISQ matrix


for (region in regionvec)  {
  
  if (region == "ALM") {
    Smat=log10(Smat_ALM+1)
    Qmat=log10(Qmat_ALM+1)
    Imat=log10(Imat_ALM+1)
    Cmat=matrix(0,length(almclid),length(almclid))
  }
  
  if (region == "VISp") {
    Smat=log10(Smat_VISp+1)
    Qmat=log10(Qmat_VISp+1)
    Imat=log10(Imat_VISp+1)
    Cmat=matrix(0,length(visclid),length(visclid))
  }
  

gmax=max(Imat,Smat,Qmat)
for (i in 1:dim(Cmat)[1]) {
  for (j in 1:dim(Cmat)[2]) {
   nval=Smat[i,j]+Qmat[i,j]+Imat[i,j]
    if (nval > 0) {
      Imat[i,j]=as.numeric((Imat[i,j]/max(Imat)))
      #Imat[i,j]=as.numeric(Imat[i,j]/gmax)
      Smat[i,j]=as.numeric((Smat[i,j]/max(Smat))) 
      #Smat[i,j]=as.numeric(Smat[i,j]/gmax)
      Qmat[i,j]=as.numeric((Qmat[i,j]/max(Qmat)))
     # Qmat[i,j]=as.numeric(Qmat[i,j]/gmax)
      Cmat[i,j]=as.character(rgb(Imat[i,j],Smat[i,j],Qmat[i,j],alpha=1))
    }
    else (Cmat[i,j]=rgb(0.75,0.75,0.75))
  }
}

fout=paste("ISQ","_",region,sep="")
write.csv(Cmat,paste("F:\\Analysis\\Neuropeptides\\ISQ\\",fout,".csv",sep=""),row.names=F)

imout=paste(fout,".tiff",sep="")
imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",imout,sep="")
#imname=paste("/volumes/PHILIPS UFD/Analysis/Neuropeptides/ISQ/",imout,sep="")

tiff(imname,width=1500,height=1500)
  if (region=="ALM") {
    Cmat_ALM=Cmat
    plot_ISQ_ALM (Cmat)
    }
  if (region=="VISp") {
    plot_ISQ_VISp (Cmat)
    Cmat_VISp=Cmat
    }
dev.off()

if (region=="ALM") {
  
  ichn="ISQ_ALM_I.tiff"
  imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",ichn,sep="")
  tiff(imname,width=1500,height=1500)
  plot_channel_ALM(Imat,"I","ALM I")
  dev.off()
  
  ichn="ISQ_ALM_S.tiff"
  imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",ichn,sep="")
  tiff(imname,width=1500,height=1500)
  plot_channel_ALM(Smat,"S","ALM S")
  dev.off()
  
  ichn="ISQ__ALM_Q.tiff"
  imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",ichn,sep="")
  tiff(imname,width=1500,height=1500)
  plot_channel_ALM(Qmat,"Q","ALM Q")
  dev.off()
  
}

if (region=="VISp") {
  
  ichn="ISQ_VISp_I.tiff"
  imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",ichn,sep="")
  tiff(imname,width=1500,height=1500)
   plot_channel_VISp(Imat,"I","VISp I")
  dev.off()
  
  ichn="ISQ_VISp_S.tiff"
  imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",ichn,sep="")
  tiff(imname,width=1500,height=1500)
  plot_channel_VISp(Smat,"S","VISp S")
  dev.off()
  
  ichn="ISQ_VISp_Q.tiff"
  imname=paste("F:\\Analysis\\Neuropeptides\\ISQ\\",ichn,sep="")
  tiff(imname,width=1500,height=1500)
  plot_channel_VISp(Qmat,"Q","VISp Q")
  dev.off()
}
 
} # end region

#################################################################################         
# ISQ matrix pooled to cell types, pool the I, S, Q matrices individually first 
#################################################################################

Cmat_pool=matrix(0,length(ctcolors),length(ctcolors))

for (region in regionvec) {
  
  if (region=="ALM") {
    
    Smat_pool=poolMat(as.matrix(Smat_ALM),region)
    Qmat_pool=poolMat(as.matrix(Qmat_ALM),region)
    Imat_pool=poolMat(as.matrix(Imat_ALM),region)
    Smat_ALM_pool=log10(Smat_pool+1)
    Qmat_ALM_pool=log10(Qmat_pool+1)
    Imat_ALM_pool=log10(Imat_pool+1)
    
  }
  
  if (region=="VISp") {
    
    Smat_pool=poolMat(as.matrix(Smat_VISp),region)
    Qmat_pool=poolMat(as.matrix(Qmat_VISp),region)
    Imat_pool=poolMat(as.matrix(Imat_VISp),region)
    Smat_VISp_pool=log10(Smat_pool+1)
    Qmat_VISp_pool=log10(Qmat_pool+1)
    Imat_VISp_pool=log10(Imat_pool+1)
  }
  
  
  
  # write out pooled matrices
  
  write.csv(Imat_ALM_pool,paste("F:\\Analysis\\Neuropeptides\\ISQ\\Imat_ALM_pool.csv",sep=""))  
  write.csv(Smat_ALM_pool,paste("F:\\Analysis\\Neuropeptides\\ISQ\\Smat_ALM_pool.csv",sep="")) 
  write.csv(Qmat_ALM_pool,paste("F:\\Analysis\\Neuropeptides\\ISQ\\Qmat_ALM_pool.csv",sep="")) 
  
  write.csv(Imat_VISp_pool,paste("F:\\Analysis\\Neuropeptides\\ISQ\\Imat_VISp_pool.csv",sep=""))
  write.csv(Smat_VISp_pool,paste("F:\\Analysis\\Neuropeptides\\ISQ\\Smat_VISp_pool.csv",sep=""))
  write.csv(Qmat_VISp_pool,paste("F:\\Analysis\\Neuropeptides\\ISQ\\Qmat_VISp_pool.csv",sep=""))
  
  if (region=="ALM")  {
    mpal = colorRampPalette(c(rgb(0,0,0),rgb(1,0,0)),bias=1)(6)
    pheatmap(Imat_ALM_pool,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
             main="ALM Pooled I",show_rownames=FALSE,show_colnames=FALSE,legend=F)
    mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,1,0)),bias=1)(6)
    pheatmap(Smat_ALM_pool,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
             main="ALM Pooled S",show_rownames=FALSE,show_colnames=FALSE,legend=F)
    mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,0,1)),bias=1)(6)
    pheatmap(Qmat_ALM_pool,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
             main="ALM Pooled Q",show_rownames=FALSE,show_colnames=FALSE,legend=F)
  }
  
  
  if (region=="VISp")  {
    mpal = colorRampPalette(c(rgb(0,0,0),rgb(1,0,0)),bias=1)(6)
    pheatmap(Imat_VISp_pool,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
             main="VISp Pooled I",show_rownames=FALSE,show_colnames=FALSE,legend=F)
    mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,1,0)),bias=1)(6)
    pheatmap(Smat_VISp_pool,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
             main="VISp Pooled S",show_rownames=FALSE,show_colnames=FALSE,legend=F)
    mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,0,1)),bias=1)(6)
    pheatmap(Qmat_VISp_pool,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
             main="VISp Pooled Q",show_rownames=FALSE,show_colnames=FALSE,legend=F)
  }
  
 
  
  for (i in 1:dim(Cmat_pool)[1]) {
    for (j in 1:dim(Cmat_pool)[2]) {
     # nval=Smat_pool[i,j]+Qmat_pool[i,j]+Imat_pool[i,j]
      if (nval > 0) {
        
        if (region == "ALM") {
          Imat_pool[i,j]=as.numeric((Imat_ALM_pool[i,j]/max(Imat_ALM_pool)))
          Smat_pool[i,j]=as.numeric((Smat_ALM_pool[i,j]/max(Smat_ALM_pool))) 
          Qmat_pool[i,j]=as.numeric((Qmat_ALM_pool[i,j]/max(Qmat_ALM_pool)))
        }
        if (region == "VISp") {
          Imat_pool[i,j]=as.numeric((Imat_VISp_pool[i,j]/max(Imat_VISp_pool)))
          Smat_pool[i,j]=as.numeric((Smat_VISp_pool[i,j]/max(Smat_VISp_pool))) 
          Qmat_pool[i,j]=as.numeric((Qmat_VISp_pool[i,j]/max(Qmat_VISp_pool)))
        }
        Cmat_pool[i,j]=as.character(rgb(Imat_pool[i,j],Smat_pool[i,j],Qmat_pool[i,j],alpha=1))
      }
      else (Cmat_pool[i,j]=rgb(0.75,0.75,0.75))
    }
  }
  
  
  # ISQ pooled visualization 
  
  line=length(ctcolors)
  col=length(ctcolors)
  
  # use many little rectangles rect(xleft, ybottom, xright, ytop)
  xleft= (0:(col - 1)/col) 
  ybottom=sort((0:(line - 1)/line) ,decreasing=T) 
  xright= (1:col/col) 
  ytop=sort((1:line/line),decreasing=T)
  
  
  par(new=F)
  plot(NULL,xlim = c(-0.2,1.2), ylim = c(-0.2,1.3),axes=F,ylab="",xlab="")
  for (i in 1:length(ctcolors)) {
    for (j in 1:length(ctcolors)) {
      rect(xleft[i],ybottom[j],xright[i],ytop[j],border = NULL , col=Cmat_pool[j,i])
    }
  }
  legend(1,1,legend=c("I","S","Q","None"),pch=15,col=c(2,3,4,rgb(0.75,0.75,0.75)),bty="n",cex=1.75)
  
  
  text(0.5,1.25,paste(region,": NPP NP-GPCR Interactions",sep=""),font=2,cex=2)
  
  # text and annotations
  
  xlt=xleft
  ybt=rep(1.0,12)
  ytt=rep(ytop[1]+2*0.08333,12)
  xrt=xright
  tposx=xlt+0.08333/2
  tposy=ybt+0.08333
  
  for (i in 1:length(ctcolors)) {
    for (j in 1:length(ctcolors)) {
      rect(xlt[i],ybt[j],xrt[i],ytt[j],border = NULL, col=ctcolors[i]) 
    }
  }
  text(tposx,tposy,ctpool,font=1,srt=90)
  
  
  xls= rep(-2*0.08333,12)
  xrs=rep(0,12)  
  yts=ytop
  ybs=ybottom
  
  tposx= -0.08333
  tposy=ybs+0.08333/2
  
  for (i in 1:length(ctcolors)) {
    for (j in 1:length(ctcolors)) {
      rect(xls[i],ybs[j],xrs[i],yts[j],border = NULL, col=ctcolors[j]) 
    }
  }
  text(tposx,tposy,ctpool,font=1,srt=0)
  
  
}
#####################    
# Plotting functions
#####################

matMax <- function(mtx)
{
  colmn <- which(mtx == max(mtx)) %/% nrow(mtx) + 1
  row <- which(mtx == max(mtx)) %% nrow(mtx)
  return( matrix(c(row, colmn), 1))
}


# helper function for pooling

poolMat<-function(inmat,reg) {
  
  poolmat=matrix(0,length(ctpool),length(ctpool),dimnames=list(ctpool,ctpool))
  
  for (j in 1:length(ctpool)) { 
    for (k in 1:length(ctpool)) {
      
      if (reg == "ALM")  {suml=almlist}
      if (reg == "VISp") {suml=vislist}
      
      poolmat[j,k] = mean(inmat[suml[[j]],suml[[k]]])
    }
  }
  
  return (poolmat)
  
}




plot_VISp<-function(datamat,titin,bval)  {
  
  
  dendid=visclust[,1]   
  clist=as.character(visclust[,4])
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(paste("CT",visdend,sep=""))
  )
  rownames(annotation_col)=paste("CT",visdend,sep="")
  
  
  names(clist)=paste("CT",visdend,sep="")
  
  ann_colors = list(CellType = clist)
  
  colnames(datamat)=as.character(annotation_col[,1])
  
  cccolors <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")),bias=bval)(10)
  
  
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=20,fontsize_col=7,fontsize_row=7,border_color=9,main=titin,annotation_col=annotation_col,
           annotation_colors = ann_colors,show_colnames=F,annotation_legend=F,color=cccolors,legend=FALSE)
  
  
  
  return()
  

}




plot_ALM<-function(datamat,titin,bval)  {
   
  
  dendid=almclust[,1]   
  clist=as.character(almclust[,4])
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(paste("CT",almdend,sep=""))
  )
  rownames(annotation_col)=paste("CT",almdend,sep="")
  
  
  names(clist)=paste("CT",almdend,sep="")
  
  ann_colors = list(CellType = clist)
  
  colnames(datamat)=as.character(annotation_col[,1])
  
  
  # For a heatmap from green to red (grid of 128 colors here)
  cccolors <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")),bias=bval)(10)
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_col=7,fontsize_row=10,border_color=9,main=titin,annotation_col=annotation_col,
           annotation_colors = ann_colors,show_colnames=F,annotation_legend=F,color=cccolors)
  
  
  
  return()
  
}



plot_All<-function(datamat,titin,bval)  {
  
  dendid=allclust[,1]   
  clist=as.character(allclust[,4])
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(paste("CT",alldend,sep=""))
  )
  rownames(annotation_col)=paste("CT",alldend,sep="")
  
  
  names(clist)=paste("CT",alldend,sep="")
  
  ann_colors = list(CellType = clist)
  
  colnames(datamat)=as.character(annotation_col[,1])
  
  
  # For a heatmap from green to red (grid of 128 colors here)
  cccolors <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")),bias=bval)(100)
  #cccolors <- colorRampPalette(c("darkblue", "dodgerblue", "gray80", "orange", "orangered"),bias=bval)(100)
    
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_col=7,fontsize_row=10,border_color=9,main=titin,annotation_col=annotation_col,
           annotation_colors = ann_colors,show_colnames=F,annotation_legend=F,color=cccolors)
  
  
  
  return()
  
}



plot_All_pooled<-function(datamat,titin,bval)  {
  
  
  poolmat=matrix(0,dim(datamat)[1],length(ctpool),dimnames=list(rownames(datamat),ctpool))
  
  for (j in 1:dim(datamat)[1]) { 
    for (k in 1:length(ctpool)) {
      poolmat[j,k] = mean(datamat[j,alllist[[k]]])
    }
  }
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(ctpool)
  )
  rownames(annotation_col)= ctpool
  
  ann_colors = list(CellType = ctcolors)
  
  # For a heatmap from green to red (grid of 128 colors here)
  cccolors <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")),bias=bval)(100)
  #cccolors <- colorRampPalette(c("darkblue", "dodgerblue", "gray80", "orange", "orangered"),bias=bval)(100)
  
  
   # imname=paste("/volumes/PHILIPS UFD/Analysis/Neuropeptides/",titin,".png",sep="")
  #  png(imname,width=800,height=500)
      pheatmap(poolmat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_col=28,fontsize_row=28,fontsize=24,border_color=9,main=titin,annotation_col=annotation_col,
           annotation_colors = ann_colors,show_colnames=T,annotation_legend=F,color=cccolors)
   # dev.off()
  
  
  return()
  
}


plot_Interaction_ALM<-function(datamat,titin, gamma)  {
  
  
  rownames(datamat)=paste("CT",almdend,sep="")
  colnames(datamat)=paste("CT",almdend,sep="")
  
  # includes modification for adding two horizontal and vertical lines for alignment with VISp
   
    dendid=almclust[,1]  
    almclid=1:84
   
    
    # now fix the annotations
  
    dendid=almclust[,1]   
    clist=as.character(almclust[,4]) 
    names(clist)=paste("CT",almdend,sep="")
    
    cvec=NULL
    for (i in 1:12) {
      cvec=c(cvec,rep(i,length(almlist[[i]])))
    }
    
    ctvec=ctcolors
    names(ctvec)=c(1:12)
    
    
    # Generate annotations for rows and columns
    annotation_col = data.frame(
      CellType = factor(paste("CT",almdend,sep="")),
      CellClass = cvec
    )
    rownames(annotation_col)=paste("CT",almdend,sep="")
    
    ann_colors = list(CellType = clist, CellClass = ctvec)
    
    clist=dendcolor[dendid]
    
    
    colnames(datamat)=as.character(annotation_col[,1])
    rownames(datamat)=as.character(annotation_col[,1])
  

    bval=1
    gamma=1
    mpal = colorRampPalette(brewer.pal(8,"PuBuGn"),bias=bval)(6)
    gammaC=1/gamma
    corcol=(col2rgb(mpal)/255)^gammaC
    mpal=rgb(t(corcol)) 
  
  
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=rev(mpal),
           main=titin,annotation_col=annotation_col,annotation_row=annotation_col,annotation_colors = ann_colors,annotation_legend=F,
           show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=F,annotation_names_col=F,legend=F)
  
} 


plot_Interaction_VISp<-function(datamat,titin,gamma)  {
  
  
  rownames(datamat)=paste("CT",visdend,sep="")
  colnames(datamat)=paste("CT",visdend,sep="")
  
  dendid=visclust[,1]   
  clist=as.character(visclust[,4]) 
  names(clist)=paste("CT",visdend,sep="")
  
  cvec=NULL
  for (i in 1:12) {
    cvec=c(cvec,rep(i,length(vislist[[i]])))
  }
  
  ctvec=ctcolors
  names(ctvec)=c(1:12)
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(paste("CT",visdend,sep="")),
    CellClass = cvec
  )
  rownames(annotation_col)=paste("CT",visdend,sep="")
  
  ann_colors = list(CellType = clist, CellClass = ctvec)
  
  colnames(datamat)=as.character(annotation_col[,1])
  rownames(datamat)=as.character(annotation_col[,1])
  
  
  bval=1
  gamma=1
  mpal = colorRampPalette(brewer.pal(8,"PuBuGn"),bias=bval)(6)
  gammaC=1/gamma
  corcol=(col2rgb(mpal)/255)^gammaC
  mpal=rgb(t(corcol)) 
  
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=rev(mpal),
           main=titin,annotation_col=annotation_col,annotation_row=annotation_col,annotation_colors = ann_colors,annotation_legend=F,
           show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=F,annotation_names_col=F,legend=F)
  
  
} 


plot_ISQ_ALM <-function (Cmat) {
  
# tricolor intensity map by I,S, Q type


line=length(almclid)
col=length(almclid)

# use many little rectangles rect(xleft, ybottom, xright, ytop)
xleft= (0:(col - 1)/col) 
ybottom=sort((0:(line - 1)/line) ,decreasing=T) 
xright= (1:col/col) 
ytop=sort((1:line/line),decreasing=T)


par(new=F)
plot(NULL,xlim = c(-0.2,1.2), ylim = c(-0.2,1.2),axes=F,ylab="",xlab="")
for (i in 1:length(almclid)) {
  for (j in 1:length(almclid)) {
    rect(xleft[i],ybottom[j],xright[i],ytop[j],border = NULL , col=Cmat[j,i])
  }
}
legend(1,1,legend=c("I","S","Q","None"),pch=15,col=c(2,3,4,rgb(0.75,0.75,0.75)),bty="n",cex=1.5)

text(0.5,1.1,"ALM: NPP NP-GPCR Interactions",font=2,cex=2)

}  


plot_ISQ_VISp <-function (Cmat) {
  
  # tricolor intensity map by I,S, Q type
  
  
  line=length(visclid)
  col=length(visclid)
  
  # use many little rectangles rect(xleft, ybottom, xright, ytop)
  xleft= (0:(col - 1)/col) 
  ybottom=sort((0:(line - 1)/line) ,decreasing=T) 
  xright= (1:col/col) 
  ytop=sort((1:line/line),decreasing=T)
  
  
  par(new=F)
  plot(NULL,xlim = c(-0.2,1.2), ylim = c(-0.2,1.2),axes=F,ylab="",xlab="")
  for (i in 1:length(visclid)) {
    for (j in 1:length(visclid)) {
      rect(xleft[i],ybottom[j],xright[i],ytop[j],border = NULL , col=Cmat[j,i])
    }
  }
  legend(1,1,legend=c("I","S","Q","None"),pch=15,col=c(2,3,4,rgb(0.75,0.75,0.75)),bty="n",cex=1.5)
  
  text(0.5,1.1,"VISp: NPP NP-GPCR Interactions",font=2,cex=2)
  
  
}  

plot_channel_ALM<-function (datamat,chnl,ititle) {
  
  # take single channel data and index and render
  
  if (chnl=="I") {mpal = colorRampPalette(c(rgb(0,0,0),rgb(1,0,0)),bias=1)(6)}
  if (chnl=="S") {mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,1,0)),bias=1)(6)}
  if (chnl=="Q") {mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,0,1)),bias=1)(6)}
  
  
  
  rownames(datamat)=paste("CT",almdend,sep="")
  colnames(datamat)=paste("CT",almdend,sep="")
  
  # includes modification for adding two horizontal and vertical lines for alignment with VISp
  
  dendid=almclust[,1]  
  almclid=1:84
  
  
  # now fix the annotations
  
  dendid=almclust[,1]   
  clist=as.character(incolors[dendid,3])
  names(clist)=paste("CT",almdend,sep="")
  
  cvec=NULL
  for (i in 1:12) {
    cvec=c(cvec,rep(i,length(almlist[[i]])))
  }
  
  ctvec=ctcolors
  names(ctvec)=c(1:12)
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(paste("CT",almdend,sep="")),
    CellClass = cvec
  )
  rownames(annotation_col)=paste("CT",almdend,sep="")
  
  ann_colors = list(CellType = clist, CellClass = ctvec)
  
  clist=dendcolor[dendid]
  
  
  colnames(datamat)=as.character(annotation_col[,1])
  rownames(datamat)=as.character(annotation_col[,1])
  
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
           main=ititle,annotation_col=annotation_col,annotation_row=annotation_col,annotation_colors = ann_colors,annotation_legend=F,
           show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=F,annotation_names_col=F,legend=F)
  

  
}




plot_channel_VISp<-function (datamat,chnl,ititle) {
  
  # take single channel data and index and render
  
  if (chnl=="I") {mpal = colorRampPalette(c(rgb(0,0,0),rgb(1,0,0)),bias=1)(6)}
  if (chnl=="S") {mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,1,0)),bias=1)(6)}
  if (chnl=="Q") {mpal = colorRampPalette(c(rgb(0,0,0),rgb(0,0,1)),bias=1)(6)}
  
 
  rownames(datamat)=paste("CT",visdend,sep="")
  colnames(datamat)=paste("CT",visdend,sep="")
  
  dendid=visclust[,1]   
  clist=as.character(incolors[dendid,3])
  names(clist)=paste("CT",visdend,sep="")
  
  cvec=NULL
  for (i in 1:12) {
    cvec=c(cvec,rep(i,length(vislist[[i]])))
  }
  
  ctvec=ctcolors
  names(ctvec)=c(1:12)
  
  
  # Generate annotations for rows and columns
  annotation_col = data.frame(
    CellType = factor(paste("CT",visdend,sep="")),
    CellClass = cvec
  )
  rownames(annotation_col)=paste("CT",visdend,sep="")
  
  ann_colors = list(CellType = clist, CellClass = ctvec)
  
  colnames(datamat)=as.character(annotation_col[,1])
  rownames(datamat)=as.character(annotation_col[,1])
   
  
  pheatmap(datamat,cluster_cols=FALSE,cluster_rows=FALSE,fontsize_row=9,fontsize_col=11,fontsize=24,border_color=9,color=mpal,
           main=ititle,annotation_col=annotation_col,annotation_row=annotation_col,annotation_colors = ann_colors,annotation_legend=F,
           show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=F,annotation_names_col=F,legend=F)
  
}

