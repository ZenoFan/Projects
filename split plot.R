Metadata=read.csv(file.choose(),check.names = F,header=T)
str(Metadata)
library(car)
library(lsmeans)
library(multcompView)
library(lme4)
library(plyr)
library(agricolae)
library(xtable)
Metadata$Genotype=as.factor(Metadata$Genotype)
Metadata$Block=as.factor(Metadata$Block)
Metadata$Time=as.factor(Metadata$Time)
Metadatanew=Metadata

leveneTest(X1.Benzylglucopyranoside ~ Genotype*Trt, data = Metadatanew)

#####loop######
options(contrasts = c("contr.sum","contr.poly"))
options(scipen=999)
twowaytable=as.matrix(data.frame(matrix(ncol = 5, nrow = 0)))
colnames(twowaytable)=c("Genotype","Trt","lsmean","SE",".group")
threewaytable=as.matrix(data.frame(matrix(ncol = 6, nrow = 0)))
colnames(threewaytable)=c("Genotype","Trt","Time","lsmean","SE",".group")
ANOVAtable=as.matrix(data.frame(matrix(ncol = 2, nrow = 0)))
colnames(ANOVAtable)=c("F value",  "Pr(>F)") 

metabnumber=0
Pvaltrt=data.frame(matrix(ncol = 2, nrow = 142))
colnames(Pvaltrt)=c("name",  "Pr")
Pvaltime=data.frame(matrix(ncol = 2, nrow = 142))
colnames(Pvaltime)=c("name",  "Pr") 
#which(colSums(is.na(Metadatanew[,5:ncol(Metadatanew)]))<5)
for (i in which(colSums(is.na(Metadatanew[,5:ncol(Metadatanew)]))<5)){
  #model=with(Metadatanew,ssp.plot(Block,Trt,Genotype,Time,Metadatanew[,i+4]))
model=aov(Metadatanew[,i+4] ~  Trt*Genotype*Time, data=Metadatanew)
  metabnumber=metabnumber+1
  model=aov(Metadatanew[,i+4] ~  Trt*Genotype*Time, data=Metadatanew)
  Pvaltrt[i,2]=round(Anova(model,test="F",type=3,singular.ok =T)[2,4],5)
  Pvaltrt[i,1]=colnames(Metadatanew[i+4])
  Pvaltime[i,2]=round(Anova(model,test="F",type=3,singular.ok =T)[4,4],5)
  Pvaltime[i,1]=colnames(Metadatanew[i+4])
  fitmodel=lmer(Metadatanew[,i+4]~Trt*Genotype*Time+(1|Block:Trt:Genotype), data=Metadatanew)
  sample.lsm1 <- lsmeans(fitmodel, pairwise~Genotype:Trt,adjust="tukey",mode="kenward-roger")
  Pvalcompt=as.data.frame(cld(sample.lsm1,alpha=.05,Letters=letters))
  #Pvalcomp=cbind(metab=colnames(Metadatanew[i+4]),Pvalcompt[,c(1:4,8)])
  Pvalcomp=rbind(c(colnames(Metadatanew[i+4]),NA,NA,NA,NA),as.matrix(Pvalcompt[,c(1:4,8)]))
  twowaytable=rbind(twowaytable,Pvalcomp)
  
  sample.lsm2 <- lsmeans(fitmodel, pairwise~Genotype:Trt:Time,adjust="tukey",mode="kenward-roger")
  fortablet=as.data.frame(cld(sample.lsm2,alpha=.05,Letters=letters))
  #fortable=cbind(metab=colnames(Metadatanew[i+4]),fortablet[,c(1:5,9)])
  fortable=rbind(c(colnames(Metadatanew[i+4]),NA,NA,NA,NA,NA),as.matrix(fortablet[,c(1:5,9)]))
  threewaytable=rbind(threewaytable,fortable)
  
  #Pvaltrtemp=cbind(metab=colnames(Metadatanew[i+4]),pvaltrt,pvalinter)
  #Pval2and3way=rbind(Pval2and3way,Pvaltrtemp)
  metabtitle=t(as.matrix(c(colnames(Metadatanew[i+4]),NA)))
  colnames(metabtitle)=c("F value",  "Pr(>F)") 
  #ANOVAtabletemp=rbind(metabtitle,as.matrix(round(model$ANOVA[-c(3,6,11),4:5],5)))
  ANOVAtabletemp=rbind(metabtitle,as.matrix(round(Anova(model,test="F",type=3,singular.ok =T)[,3:4],5)))
  ANOVAtable=rbind(ANOVAtable,ANOVAtabletemp)
}
ANOVAtable[is.na(ANOVAtable)]=""
twowaytable[is.na(twowaytable)]=" "
threewaytable[is.na(threewaytable)]=" "
write.csv(ANOVAtable, "RootANOVAtable.csv")
write.csv(twowaytable, "Roottwowaytable.csv", row.names=FALSE)
write.csv(threewaytable, "Rootthreewaytable.csv", row.names=FALSE)

selectPtrt=as.data.frame(Pvaltrt[order(Pvaltrt$Pr),])[1:30,1]
selectPtime=as.data.frame(Pvaltime[order(Pvaltime$Pr),])[1:30,1]
selectMetadatanew=Metadatanew[,(colnames(Metadatanew)%in%selectPtrt)|(colnames(Metadatanew)%in%selectPtime)]
selectMetadatanew=cbind(Metadatanew[,1:4],selectMetadatanew)

selecttwowaytable=as.matrix(data.frame(matrix(ncol = 5, nrow = 0)))
colnames(selecttwowaytable)=c("Genotype","Trt","lsmean","SE",".group")
selectthreewaytable=as.matrix(data.frame(matrix(ncol = 6, nrow = 0)))
colnames(selectthreewaytable)=c("Genotype","Trt","Time","lsmean","SE",".group")
metabnumber1=0

for (i in (1:ncol(selectMetadatanew))){
  fitmodel=lmer(selectMetadatanew[,i+4]~Trt*Genotype*Time+(1|Block:Trt:Genotype), data=selectMetadatanew)
  sample.lsm1 <- lsmeans(fitmodel, pairwise~Genotype:Trt,adjust="tukey",mode="kenward-roger")
  Pvalcompt=as.data.frame(cld(sample.lsm1,alpha=.05,Letters=letters))
  Pvalcomp=rbind(c(colnames(selectMetadatanew[i+4]),NA,NA,NA,NA),as.matrix(Pvalcompt[,c(1:4,8)]))
  selecttwowaytable=rbind(selecttwowaytable,Pvalcomp)
  
  metabnumber1=metabnumber1+1
  
  sample.lsm2 <- lsmeans(fitmodel, pairwise~Genotype:Trt:Time,adjust="tukey",mode="kenward-roger")
  fortablet=as.data.frame(cld(sample.lsm2,alpha=.05,Letters=letters))
   fortable=rbind(c(colnames(selectMetadatanew[i+4]),NA,NA,NA,NA,NA),as.matrix(fortablet[,c(1:5,9)]))
  selectthreewaytable=rbind(selectthreewaytable,fortable)
}
selecttwowaytable[is.na(selecttwowaytable)]=" "
selectthreewaytable[is.na(selectthreewaytable)]=" "
write.csv(selecttwowaytable, "selectroottwowaytable.csv", row.names=FALSE)
write.csv(selectthreewaytable, "selectrootthreewaytable.csv", row.names=FALSE)
