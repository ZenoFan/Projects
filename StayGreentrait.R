library(sommer)
library(lmerTest)
library(pbkrtest)
library(car)
library(lsmeans)
options(scipen=999)
#####4 sampling dates centered on ####################
Mdata = read.csv(file.choose(), header=TRUE, sep=",",stringsAsFactors=FALSE,check.names = F)
GDDdata = read.csv(file.choose(), header=TRUE, sep=",",stringsAsFactors=FALSE,check.names = F)
colnames(GDDdata)[4]="Days"
GDDdataCitra=GDDdata[GDDdata[,1]=='Citra',]
rownames(GDDdataCitra)=NULL
GDDdataQuincy=GDDdata[GDDdata[,1]=='Quincy',]
rownames(GDDdataQuincy)=NULL


lookupCitra=c(73,86,96,107,114,122,128,136)
lookupQuincy=c(76,87,94,110,118,125,131)

lookupCitraFORGDD=lookupCitra+18
lookupQuincyFORGDD=lookupQuincy+18
Mdata$NDVI_1=NA
Mdata$NDVI_2=NA
Mdata$NDVI_3=NA
Mdata$NDVI_4=NA
Mdata$NDVI_5=NA
Mdata$NDVI_6=NA

Mdata$GDD_1=NA
Mdata$GDD_2=NA
Mdata$GDD_3=NA
Mdata$GDD_4=NA
Mdata$GDD_5=NA
Mdata$GDD_6=NA
Mdata$Miday=Mdata$`Days to heading`+(Mdata$`Days to Maturity`-Mdata$`Days to heading`)/2
Mdata=Mdata[,c(1:7,40,8:39)]
Mdata$Miday=Mdata$Miday+18

for (i in 1:nrow(Mdata)){
  if(is.na(Mdata[i,7])|is.na(Mdata[i,18])){next}else{
if(Mdata[i,colnames(Mdata)=='Location']=='Citra')
{
  k=abs(lookupCitraFORGDD-Mdata[i,8])
  l=(k==max(k[k!=max(k)])|k==max(k))
  m=lookupCitraFORGDD[!l]
  if(length(m)==5) o=c(T,F,F,F,F,F,F,T) else o=l
  Mdata[i,29:34]=Mdata[i,18:25][!o]
  if(length(m)==5) n=lookupCitraFORGDD[!o] else n=m
  Mdata[i,35:40]= GDDdataCitra[n,5]
}else 
{
  k=abs(lookupQuincyFORGDD-Mdata[i,8])
  l=(k==max(k))
  m=lookupQuincyFORGDD[!l]
  if(length(m)==5) o=ifelse(k[2]>k[6],c(T,F,F,F,F,F,F),c(F,F,F,F,F,F,T)) else o=l
  Mdata[i,29:34]=Mdata[i,18:24][!o]
  if(length(m)==5) n=lookupQuincyFORGDD[!o] else n=m
  Mdata[i,35:40]= GDDdataCitra[n,5]
}
}
  print(i)
}

#Mdata[which(Mdata[,32]==Mdata[,28]),32]=Mdata[which(Mdata[,32]==Mdata[,28]),24]

Mdata$DTMnew=Mdata$`Days to Maturity`+18
Mdata$MTT=NA
for (i in 1:nrow(Mdata)){
  if(Mdata[i,colnames(Mdata)=='Location']=='Citra')
  {
  Mdata$MTT[i]=GDDdataCitra[GDDdataCitra[,4]==Mdata$DTMnew[i],5]  
  }else 
  {
  Mdata$MTT[i]=GDDdataQuincy[GDDdataQuincy[,4]==Mdata$DTMnew[i],5]  
  }
  print(i)
}

Mdata$NDVI_1=(Mdata$NDVI_1)/100
Mdata$NDVI_2=(Mdata$NDVI_2)/100
Mdata$NDVI_3=(Mdata$NDVI_3)/100
Mdata$NDVI_4=(Mdata$NDVI_4)/100
Mdata$NDVI_5=(Mdata$NDVI_5)/100
Mdata$NDVI_6=(Mdata$NDVI_6)/100
#Mdata = read.csv(file.choose(), header=TRUE, sep=",",stringsAsFactors=FALSE,check.names = F)
Mdata$RS=NA
Mdata$Stg=NA
for (i in 1:nrow(Mdata)){
  if (is.na(Mdata[i,7])|is.na(Mdata[i,18])){next}else{
    y=t(Mdata[i,29:34])
    
    x=t(Mdata[i,35:40])
    lm1 = lm(y ~ as.vector(x))
    RS=summary(lm1)$coefficients[2]
    Stg=predict(lm1, data.frame(x=c(Mdata[i,42])), level = 0.95, interval = "confidence")[1]
    Mdata$RS[i]=RS
    Mdata$Stg[i]=Stg
    print(i)
  }
}

P1_1_NDVI=t(Mdata[Mdata[,1]=='NC06-19896 (Parent-2)',29:34])[,1]
P1_2_NDVI=t(Mdata[Mdata[,1]=='NC06-19896 (Parent-2)',29:34])[,2]
P1_1_TT=t(Mdata[Mdata[,1]=='NC06-19896 (Parent-2)',35:40])[,1]
P1_2_TT=t(Mdata[Mdata[,1]=='NC06-19896 (Parent-2)',35:40])[,2]
P1data=as.data.frame(cbind(c(P1_1_NDVI,P1_2_NDVI),c(P1_1_TT,P1_2_TT)))
P1data$genotype[1:6]=c("NC_Citra") 
P1data$genotype[7:12]=c("NC_Quincy")


P2_1_NDVI=t(Mdata[Mdata[,1]=='AGS2000 (Parent-1)',29:34])[,1]
P2_2_NDVI=t(Mdata[Mdata[,1]=='AGS2000 (Parent-1)',29:34])[,2]
P2_1_TT=t(Mdata[Mdata[,1]=='AGS2000 (Parent-1)',35:40])[,1]
P2_2_TT=t(Mdata[Mdata[,1]=='AGS2000 (Parent-1)',35:40])[,2]
P2data=as.data.frame(cbind(c(P2_1_NDVI,P2_2_NDVI),c(P2_1_TT,P2_2_TT)))
P2data$genotype[1:6]=c("AGS_Citra") 
P2data$genotype[7:12]=c("AGS_Quincy")

Pdata=rbind(P1data,P2data)
rownames(Pdata)=NULL

library(ggplot2)
g = ggplot(Pdata, aes(x=V2, y=V1, col=genotype))+geom_point()+
  #geom_smooth(se = FALSE,method = "lm")
geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)
print(g)


plot(P1data[,1],P1data[,3],col="red")
po(P1data[,2],P1data[,4],col="red")
abline(lm(P1y1~P1x1,data=P1data))

plot(Mdata$Yield,Mdata$Stg)
plot(Mdata$Yield,Mdata$RS)

write.csv(Data_4_16,"Muhsin_4_16.csv")
Data_4_16=Newdata2[-c(194,282,251,125),]
#####Model ####################
library(sommer)
library(lmerTest)
library(pbkrtest)
library(car)
#Malldata <- read.csv(file.choose(), header=TRUE, sep=",",stringsAsFactors=FALSE,check.names = F)
#as.data.frame(mygeno)[1:4]
#write.csv(rownames(as.data.frame(mygeno)),"dipgenonames.csv")
#mygeno=as.data.frame(mygeno)
#Dipgeno=mygeno[rownames(mygeno)%in%Dipdata[,9],]
#Dippheno=Dipdata[Dipdata[,9]%in%rownames(Dipgeno),]
Mdata$row=as.integer(Mdata$row)
Mdata$column=as.integer(Mdata$column)
Mdata$rowf=as.factor(Mdata$row)
Mdata$Columnf=as.factor(Mdata$column)

Mdata$Block=as.factor(Mdata$Block)
Mdata$newBlock=paste(Mdata$Location,Mdata$Block,sep ="_")
Mdata$newBlock=as.factor(Mdata$newBlock)
Mdata$Location=as.factor(Mdata$Location)
Mdata$Entry=as.factor(Mdata$Entry)
#str(MLRdata)
#Dipphenonew <- aggregate(Dippheno[,3:ncol(Dippheno)], by = list(PEDIGREE = Dippheno$PEDIGREE), FUN = mean)
#rownames(Dipphenonew)=Dipphenonew[,1]
#### create the variance-covariance matrix
#A <- A.mat(Dipgeno) # additive relationship matrix
#### look at the data and fit the model
head(Mdata)
mix1 <- mmer2(Pred.Yield~1,
              random=~Entry
              +Location+Entry:Location
              + rowf + Columnf
              +spl2D(row,column, at=Location),
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Newdata2[!Newdata2[,1]=="SS8641",])
mix2 <- mmer2(Pred.Stg~1,
              random=~Entry
              +Location+Entry:Location
              + rowf + Columnf
              +spl2D(row,column, at=Location),
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Newdata2[!Newdata2[,1]=="SS8641",])
mix3 <- mmer2(Pred.Yield~1,
              random=~Entry
              +Location+Entry:Location
              +newBlock,
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Newdata2[!Newdata2[,1]=="SS8641",])
fmtest <- lmer(Yield~(1|Entry)+(1|Location)+(1|Location:Entry)+(1|newBlock),data=Newdata2[!Newdata2[,1]=="SS8641",])
plot(Newdata2$Yield, Newdata2$Pred.Yield)
abline(1,1)
summary()
summary(mix1)
summary(mix3)
summary(fmtest)

fm1 <- lm(Yield~vern+Entry:Location,data=Mdata)
fm1grid <- ref.grid(fm1)
lsmeansYield=summary(fm1grid)[-1]
Newdata=merge(Mdata, lsmeansYield[,1:3], by = c("Entry","Location"),all=TRUE)
fm2 <- lm(Stg~vern+Entry:Location,data=Mdata)
fm2grid <- ref.grid(fm2)
lsmeansStg=summary(fm2grid)[-1]
Newdata1=merge(Newdata, lsmeansStg[,1:3], by = c("Entry","Location"),all=TRUE)
fm3 <- lm(RS~vern+Entry:Location,data=Mdata)
fm3grid <- ref.grid(fm3)
lsmeansRS=summary(fm3grid)[-1]
Newdata2=merge(Newdata1, lsmeansRS[,1:3], by = c("Entry","Location"),all=TRUE)
colnames(Newdata2)[48:50]=c("Pred.Yield","Pred.Stg","Pred.RS")

plot(Newdata2$Pred.RS,Newdata2$Pred.Yield)
plot(Newdata2$RS,Newdata2$Yield)
abline(lm(Pred.Yield~ Pred.RS,data=Newdata2))
summary(lm(Pred.Yield~ Pred.RS,data=Newdata2))

plot(Newdata2$Pred.Stg,Newdata2$Pred.Yield)
plot(Newdata2$Stg,Newdata2$Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Newdata2))
summary(lm(Pred.Yield~ Pred.Stg,data=Newdata2))
summary(lm(Yield~ Stg,data=Newdata2))

g = ggplot(Newdata2, aes(x=Pred.Stg, y=Pred.Yield, col=Location))+geom_point()+
  geom_smooth(se = FALSE,method = "lm")
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)
print(g)

plot(Newdata2[Newdata2[,2]=='Quincy',]$Pred.Stg,Newdata2[Newdata2[,2]=='Quincy',]$Pred.Yield)
plot(Newdata2[Newdata2[,2]=='Quincy',]$Stg,Newdata2[Newdata2[,2]=='Quincy',]$Yield)
abline(lm(Newdata2[Newdata2[,2]=='Quincy',]$Pred.Yield~ Newdata2[Newdata2[,2]=='Quincy',]$Pred.Stg,data=Newdata2))
summary(lm(Newdata2[Newdata2[,2]=='Quincy',]$Pred.Yield~ Newdata2[Newdata2[,2]=='Quincy',]$Pred.Stg,data=Newdata2))
summary(lm(Yield~ Stg,data=Newdata2[Newdata2[,2]=='Quincy',]))

plot(Newdata2[Newdata2[,2]=='Citra',]$vern,Newdata2[Newdata2[,2]=='Citra',]$Yield)

summary(mix1)
YieldBLUP=as.data.frame(mix1$u.hat$Entry)
YieldBLUP$Entry=rownames(YieldBLUP)
StgBLUP=as.data.frame(mix2$u.hat$Entry)
StgBLUP$Entry=rownames(StgBLUP)
Mcorrdata=merge(YieldBLUP,StgBLUP,by=2)

plot(Mcorrdata$Pred.Stg,Mcorrdata$Pred.Yield)
plot(Mcorrdata$Stg,Mcorrdata$prediction)
write.csv(Newdata2,"Mushin_pred_4_13.csv")
#### Quincy
Mquincy=Mdata[Mdata[,6]=="Quincy",]
mix1 <- mmer2(Pred.Yield~1,
              random=~Entry
              + rowf + Columnf
              +spl2D(row,column),
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Newquincydata2)
mix2 <- mmer2(Pred.Stg~1,
              random=~Entry
              + rowf + Columnf
              +spl2D(row,column),
              rcov=~units, 
              #G=list(id=A),
              silent=T,
              data=Newquincydata2)
fm1 <- lm(Yield~vern+Entry,data=Mquincy)
fm1grid <- ref.grid(fm1)
lsmeansYield=summary(fm1grid)[-1]
Newquincydata=merge(Mquincy, lsmeansYield[,1:2], by = c("Entry"),all=TRUE)
fm2 <- lm(Stg~vern+Entry,data=Mquincy)
fm2grid <- ref.grid(fm2)
lsmeansStg=summary(fm2grid)[-1]
Newquincydata1=merge(Newquincydata, lsmeansStg[,1:2], by = c("Entry"),all=TRUE)
fm3 <- lm(RS~vern+Entry,data=Mquincy)
fm3grid <- ref.grid(fm3)
lsmeansRS=summary(fm3grid)[-1]
Newquincydata2=merge(Newquincydata1, lsmeansRS[,1:2], by = c("Entry"),all=TRUE)
colnames(Newquincydata2)[43:45]=c("Pred.Yield","Pred.Stg","Pred.RS")

plot(Newquincydata2$Pred.RS,Newquincydata2$Pred.Yield)
plot(Newquincydata2$RS,Newquincydata2$Yield)
abline(lm(Pred.Yield~ Pred.RS,data=Newquincydata2))
summary(lm(Pred.Yield~ Pred.RS,data=Newquincydata2))

plot(Newquincydata2$Pred.Stg,Newquincydata2$Pred.Yield)
plot(Newquincydata2$Stg,Newquincydata2$Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Newquincydata2))
summary(lm(Pred.Yield~ Pred.Stg,data=Newquincydata2))


plot(Mcitra$RS,Mcitra$Yield)
plot(Newdata2$Stg,Newdata2$Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Newdata2))
summary(lm(Pred.Yield~ Pred.RS,data=Newquincydata2))

plot(Mcitra$Stg,Mcitra$Yield)
plot(Newdata2$Stg,Newdata2$Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Newdata2))

summary(mix1)
YieldBLUP=as.data.frame(mix1$u.hat$Entry)
YieldBLUP$Entry=rownames(YieldBLUP)
StgBLUP=as.data.frame(mix2$u.hat$Entry)
StgBLUP$Entry=rownames(StgBLUP)
Mcorrdata=merge(YieldBLUP,StgBLUP,by=2)

plot(Mcorrdata$Pred.Stg,Mcorrdata$Pred.Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Mcorrdata))
plot(Mcorrdata$Stg,Mcorrdata$prediction)
#### Citra
Mcitra=Mdata[Mdata[,6]=="Citra",]
mix1 <- mmer2(Pred.Yield~1,
              random=~Entry
              + rowf + Columnf
              +spl2D(row,column),
              rcov=~units, 
              #G=list(id=A),
              silent=TRUE,
              data=Newcitradata2)
mix2 <- mmer2(Pred.Stg~1,
              random=~Entry
              + rowf + Columnf
              +spl2D(row,column),
              rcov=~units, 
              #G=list(id=A),
              silent=T,
              data=Newcitradata2)
fm1 <- lm(Yield~vern+Entry,data=Mcitra)
fm1grid <- ref.grid(fm1)
lsmeansYield=summary(fm1grid)[-1]
Newcitradata=merge(Mcitra, lsmeansYield[,1:2], by = c("Entry"),all=TRUE)
fm2 <- lm(Stg~vern+Entry,data=Mcitra)
fm2grid <- ref.grid(fm2)
lsmeansStg=summary(fm2grid)[-1]
Newcitradata1=merge(Newcitradata, lsmeansStg[,1:2], by = c("Entry"),all=TRUE)
fm3 <- lm(RS~vern+Entry,data=Mcitra)
fm3grid <- ref.grid(fm3)
lsmeansRS=summary(fm3grid)[-1]
Newcitradata2=merge(Newcitradata1, lsmeansRS[,1:2], by = c("Entry"),all=TRUE)
colnames(Newcitradata2)[37:39]=c("Pred.Yield","Pred.Stg","Pred.RS")

plot(Newcitradata2$Pred.RS,Newcitradata2$Pred.Yield)
plot(Newcitradata2$RS,Newcitradata2$Yield)
abline(lm(Pred.Yield~ Pred.RS,data=Newcitradata2))
summary(lm(Pred.Yield~ Pred.RS,data=Newcitradata2))

plot(Newcitradata2$Pred.Stg,Newcitradata2$Pred.Yield)
plot(Newcitradata2$Stg,Newcitradata2$Yield)
abline(lm(Yield~ Stg,data=Newcitradata2))
summary(lm(Yield~ Stg,data=Newcitradata2))

plot(Mcitra$RS,Mcitra$Yield)
plot(Newdata2$Stg,Newdata2$Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Newdata2))
summary(lm(Pred.Yield~ Pred.RS,data=Newcitradata2))

plot(Mcitra$Stg,Mcitra$Yield)
plot(Newdata2$Stg,Newdata2$Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Newdata2))

summary(mix1)
YieldBLUP=as.data.frame(mix1$u.hat$Entry)
YieldBLUP$Entry=rownames(YieldBLUP)
StgBLUP=as.data.frame(mix2$u.hat$Entry)
StgBLUP$Entry=rownames(StgBLUP)
Mcorrdata=merge(YieldBLUP,StgBLUP,by=2)

plot(Mcorrdata$Pred.Stg,Mcorrdata$Pred.Yield)
abline(lm(Pred.Yield~ Pred.Stg,data=Mcorrdata))
plot(Mcorrdata$Stg,Mcorrdata$prediction)


summary(mix1)
summary(mix2)
summary(fm1)
length(levels(Mdata$Location))
pin(mix3, h2 ~ V1 / ( V1 + (V3/2) + (V5/(2)) ) )
mix3$var.comp
  
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




