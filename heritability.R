options(contrasts=c("contr.sum","contr.poly"))
options(scipen=999)
Islam_msu=read.table("Phenotype_All LY_MSU.csv", header=T, sep=",")
str(Islam_msu)
Islam_msu$year<- factor(Islam_msu$year)
Islam_msu$rep <- factor(Islam_msu$rep)
markerFULL<-read.table("Genotypes-Imp.txt",h=TRUE) 

library(bigmemory)
library(biganalytics)
library(sommer)
######################broad sense h2##############################################
str(Islam_msu)
Islam_msu$block=paste(Islam_msu$year,sep = "_",Islam_msu$rep)
IslamSFC=as.matrix(as.data.frame(Islam_msu$SFC))

ans1=mmer2(IslamSFC~1, random=~RILno+year+RILno:year+block,data=Islam_msu,silent=TRUE)
vc=ans1$var.comp
V_G=vc[[1]]
V_GE=vc[[3]]
V_E=vc[[2]]
Ve=vc[[5]]
str(Islam_msu)
n.env <- length(levels(Islam_msu$year))
pin(ans1,h2broad~V_G/(V_G + (V_GE/n.env) + (Ve/(3*n.env))))  

######################narrow sense h2##############################################
X<-as.matrix(as.data.frame(markerFULL))
X<-t(X)
IslamY=read.table("TempAdjMeans.csv",h=TRUE,sep=',') 
colnames(IslamY)[2]="id"
#Xconvert <- apply(X,2,function(x){x[which(x == 0)] <- -1; return(x)})

library("scales")
Xconvert=X
Xconvert=rescale(Xconvert, to=c(-1,1))

#table(X!=1&X!=2&X!=0)
#Xconvert=X
#Xconvert[Xconvert==0]=-1
#Xconvert[Xconvert==1]=0
#Xconvert[Xconvert==2]=1
#Xconvert[(Xconvert!=0)&(Xconvert!=-1)&(Xconvert!=1)]=NA
#write.csv(t(Xconvert),'Xconvert.csv')


#source("https://bioconductor.org/biocLite.R")
#biocLite("snpStats")
#library(snpStats)
#showClass("ImputationRules")

#Xconvert[is.na(Xconvert)]=0.5

#sum(is.na(Xconvert))
#table(Xconvert)
rownames(Xconvert)=IslamY$id
A1=A.mat(Xconvert)
D1=D.mat(Xconvert)
E1=E.mat(Xconvert)
rownames(IslamY)=IslamY$RIL

IslamY$idd=IslamY$id
IslamY$ide=IslamY$id
ytest=as.matrix(IslamY$ELO_2010)
rownames(ytest)=IslamY$id

ans.ADE <- mmer2(ytest~1,
                 random=~g(id)+g(idd)+g(ide),
                 rcov=~units,
                 G=list(id=A1,idd=D1,ide=E1),
                 silent = TRUE, data=IslamY)
ans.A <- mmer2(ytest~1,
                 random=~g(id),
                 rcov=~units,
                 G=list(id=A1),
                 silent = TRUE, data=IslamY)

(suma <- summary(ans.ADE)$var.comp.table)
(suma <- summary(ans.A)$var.comp.table)

(h2 <- sum(suma[1,1])/sum(suma[,1]))

(H2 <- sum(suma[1:3,1])/sum(suma[,1]))
Vary<-var(ytest,na.rm=TRUE)
h2=suma[1,1]/Vary

#### calculate heritability
pin(ans.A, h1 ~ V1/(V1+V2) )
#### calculate additive variance
pin(ans.A, h1 ~ V1 )
#### calculate phenotypic variance
pin(ans.A, h1 ~ V1+V2 )
## not the same than var(y)
var(CPpheno$Yield, na.rm=TRUE)
str(ans.A)



######################narrow sense h2 with hetero##############################################
colnames(Islam_msu)[1]="id"
Islam_test=Islam_msu[Islam_msu[,1]%in%IslamY[,2],]

Islam_test$idd=Islam_test$id
Islam_test$ide=Islam_test$id
rownames(Islam_test)=NULL

ytest=as.matrix(Islam_test$ELO)
rownames(ytest)=Islam_test$id


ans.ADEnew <- mmer2(ytest~1,
                 random=~g(id)+g(idd)+g(ide),
                 rcov=~units,
                 G=list(id=A1,idd=D1,ide=E1),
                 silent = TRUE, data=Islam_test)


ans.Anew <- mmer2(ytest~1,
               random=~g(id),
               rcov=~units,
               G=list(id=A1),
               silent = TRUE, data=IslamY)

(suma <- summary(ans.ADE)$var.comp.table)
(suma <- summary(ans.A)$var.comp.table)

(h2 <- sum(suma[1,1])/sum(suma[,1]))

(H2 <- sum(suma[1:3,1])/sum(suma[,1]))
Vary<-var(ytest,na.rm=TRUE)
h2=suma[1,1]/Vary





Z1 <- model.matrix(~id-1, IslamY); colnames(Z1) <- gsub("id","",colnames(Z1))

ETA <- list(id=list(Z=Z1,K=A1),idd=list(Z=Z1,K=D1),ide=list(Z=Z1,K=E1))
ans.ADE <- mmer(Y=ytest, Z=ETA, silent = TRUE)
ans.ADE$var.comp


data(CPdata)
CPpheno$idd <-CPpheno$id; CPpheno$ide <-CPpheno$id
### look at the data
head(CPpheno)
CPgeno[1:5,1:4]
A <- A.mat(CPgeno) # additive relationship matrix
D <- D.mat(CPgeno) # dominance relationship matrix
E <- E.mat(CPgeno) # epistatic relationship matrix
ans.ADE <- mmer2(color~1,
                 random=~g(id) + g(idd) + g(ide),
                 rcov=~units,
                 G=list(id=A,idd=D,ide=E),
                 silent = TRUE, data=CPpheno)
suma <- summary(ans.ADE)$var.comp.table
(h2 <- sum(suma[1,1])/sum(suma[,1]))



library(heritability)
marker_h2(ytest, X, covariates = NULL, K, alpha = 0.05,
          eps = 1e-06, max.iter = 100, fix.h2 = FALSE, h2 = 0.5)

