library(sommer)
library(lmerTest)
library(pbkrtest)
library(car)
Dipdata <- read.csv(file.choose(), header=TRUE, sep=",",stringsAsFactors=FALSE,check.names = F)
#as.data.frame(mygeno)[1:4]

#write.csv(rownames(as.data.frame(mygeno)),"dipgenonames.csv")
mygeno=as.data.frame(mygeno)
Dipgeno=mygeno[rownames(mygeno)%in%Dipdata[,9],]
Dippheno=Dipdata[Dipdata[,9]%in%rownames(Dipgeno),]
Dippheno$row=as.integer(Dippheno$row)
Dippheno$Column=as.integer(Dippheno$Column)
Dippheno$rowf=as.factor(Dippheno$row)
Dippheno$Columnf=as.factor(Dippheno$Column)

Dippheno$Block=as.factor(Dippheno$Block)
Dippheno$newBlock=paste(Dippheno$Location,Dippheno$Block,sep ="_")
Dippheno$newBlock=as.factor(Dippheno$newBlock)
Dippheno$Location=as.factor(Dippheno$Location)
Dippheno$Pedigree=as.factor(Dippheno$Pedigree)
Dippheno$Entry=as.factor(Dippheno$Entry)
str(Dippheno)
#Dipphenonew <- aggregate(Dippheno[,3:ncol(Dippheno)], by = list(PEDIGREE = Dippheno$PEDIGREE), FUN = mean)
#rownames(Dipphenonew)=Dipphenonew[,1]
#### create the variance-covariance matrix
A <- A.mat(Dipgeno) # additive relationship matrix
#### look at the data and fit the model
head(Dippheno)
DipCitra=Dippheno[Dippheno[,1]=="Citra",]
mix1 <- mmer2(YIELD~1+Location,
              random=~Entry
              +Location+Entry:Location
              + rowf + Columnf
              +spl2D(row,Column, at=Location),
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Dippheno)
mix2 <- mmer2(YIELD~1,
              random=~Entry
              +Location+Entry:Location+newBlock,
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Dippheno)
fm1 <- lmer(YIELD~(1|Entry)+ (1|Location) + (1|newBlock)+ (1|Entry:Location),data=Dippheno)
summary(mix1)
summary(mix2)
summary(fm1)
length(levels(Dippheno$Location))
pin(mix2, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )
mix2$

#### get the spatial plots
fittedvals <- spatPlots(mix1,row = "row", range = "Column")
#### extract blups
mix1$u.hat$Entry
pin(mix1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

V1=mix1$var.comp$Pedigree
V2=mix1$var.comp$Location
V3=mix1$var.comp$units

#### genetic correlation
gvc <- mix1$var.comp$Pedigree
sd.gvc <- as.matrix(sqrt(diag(gvc)))
prod.sd <- sd.gvc %*% t(sd.gvc)
(gen.cor <- gvc/prod.sd)
#### heritability

