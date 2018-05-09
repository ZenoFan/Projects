library(ggplot2)
library(rgl)
library(pegas)
library(adegenet)
library(DMwR)
library(scrime)
library(rgl)
library(dendextend)
library(beepr)
source("https://bioconductor.org/biocLite.R")
biocLite("trio")
library(trio)

hapMap2genlight <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)]
  samples <- names(hapmap)[-1]
  loci <- row.names(hapmap)
  
  # set up conversion table
  s <- as.integer(c(0,1,2,NA))
  ac <- s
  ag <- s
  at <- s
  cg <- s
  ct <- s
  gt <- s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G")
  
  # Pull out and convert genotypes
  S <- length(samples)
  SBlist <- vector(mode="list",S)   # set up list of SNPbin objects
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=TRUE, USE.NAMES=FALSE)
    # create SNPbin object for this individual
    SBlist[[i]] <- new("SNPbin", mygen)
  }
  
  # make genlight object
  x <- new("genlight", SBlist)
  locNames(x) <- loci
  indNames(x) <- samples
  
  return(x)
}
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # Add plane
  if(show.plane) 
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
}
lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
ECOcleaned2A=read.table(file.choose(),header=T)
Famlist=as.data.frame(ECOcleaned2A[,2])
mylist=read.csv(file.choose())

#inputPed=read.pedfile(file.choose(),coded = "ATCG")

ECOcleaned2A[,2]=gsub('-', '_', ECOcleaned2A[,2],fixed=T)
ECOcleaned2A[,2]=gsub('.', '_', ECOcleaned2A[,2],fixed=T)
ECOcleaned2A=ECOcleaned2A[ECOcleaned2A[,2]%in%mylist[,2],]
DAPClist=mylist[mylist[,2]%in%ECOcleaned2A[,2],]
DAPClist=mylist
ECOcleaned2A[,2]=mylist[,2]

ECOcleaned2Anew=cbind(ECOcleaned2A[,1:6],afterHWEtest)
write.csv(ECOcleaned2Anew,file ="ECOcleaned2Anew.csv")
#########quality control#########################################################################################################################################################
ECOdata=ECOcleaned2A[,7:ncol(ECOcleaned2A)]
miss_ind=c()
for (i in 1:nrow(ECOdata)){
  miss_ind[i] = sum(is.na(ECOdata[i,]))/ncol(ECOdata)
}
hist(miss_ind)
mean(miss_ind)
newECOdata=ECOdata[miss_ind<0.2,]
miss_snp=c()
for (i in 1:ncol(ECOdata)){
  miss_snp[i] = sum(is.na(ECOdata[,i]))/ncol(ECOdata)
}
hist(miss_snp)
newECOdata=newECOdata[,miss_snp<1]
mean(miss_snp)

chisq_pval = c() 
pb <- txtProgressBar(min = 5, max = ncol(newECOdata), style = 3)
for (i in 1:ncol(newECOdata)) {
  #print(i)
  #Calculate the genotype frequencies for each of the three possible genotypes.  Note that "na.rm = T" must be used to exclude NAs
  homo_ref_sum = sum(newECOdata[,i] == 0, na.rm = T)
  hetero_sum = sum(newECOdata[,i] == 1, na.rm = T)
  homo_alt_sum = sum(newECOdata[,i] == 2,na.rm=T)
  #Get the sample size
  sample_size = sum(homo_ref_sum, hetero_sum, homo_alt_sum)
  #Calculate the allele frequencies for the two alleles.
  ref_freq = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
  alt_freq = ((homo_alt_sum * 2) + hetero_sum)/(sample_size * 2)
  #Calculate the expected genotype frequencies under Hardy-Weinberg Equilibrium.
  expected_homo_ref = ref_freq^2
  expected_hetero = 2*ref_freq*alt_freq
  expected_homo_alt = alt_freq^2
  #Perform a chisq.test to test for a difference between the expected and observed genotype frequencies.
  expected_geno_freqs = c(expected_homo_ref, expected_hetero, expected_homo_alt)
  observed_geno_freqs = c(homo_ref_sum, hetero_sum, homo_alt_sum)
  chisq_test = chisq.test(observed_geno_freqs, p = expected_geno_freqs)
  #the function chisq.test() returns a "list" of results. See help(chisq.test) under the heading "Value" for details of all results calculated by this function.
  #A "list" is a special type of object in R. In the command below, the "$" symbol is used to extract a p-value from a list of chisq.test results.
  chisq_pval[i] = chisq_test$p.value
  setTxtProgressBar(pb, i+5)
}
hist (chisq_pval)
hist(chisq_pval,breaks=200)

afterHWEtest=newECOdata[,chisq_pval>=0.05]
newmylist=mylist[rownames(afterHWEtest),]

table(afterHWEtest[rownames(afterHWEtest)==67,]==afterHWEtest[rownames(afterHWEtest)==87,])
dim(afterHWEtest)

#matching alleles
matchingdata=afterHWEtest

matchingdata=cbind(as.numeric(newmylist[,1]),afterHWEtest)
colnames(matchingdata)[1]='ID'
rownames(matchingdata)=newmylist[,4]


NaCountingforrow=matrix(NA,nrow(matchingdata),3)
for (i in 1:nrow(matchingdata)){
  NaCountingforrow[i,1]=rownames(matchingdata)[i]
  NaCountingforrow[i,2]=sum(is.na(matchingdata[i,]))
  for (j in i:nrow(matchingdata)){
    if(rownames(matchingdata)[j]==rownames(matchingdata)[i])
      NaCountingforrow[i,3]=sum(as.numeric(matchingdata[i,-1]==matchingdata[j,-1]),na.rm = T)/length(which(!is.na(matchingdata[j,-1])))
  }
  print(i)

  sum(as.numeric(matchingdata[134,-1]==matchingdata[154,-1]),na.rm = T)/length(which(!is.na(matchingdata[154,-1])))
  #Allele frequency test
  frequency=c()
  for (i in 1:ncol(mydata2)){
    frequency[i]=round(sum(mydata2[,i])/90,2)  
    print(i)
  }
  newmydata2=mydata2[,frequency>0.2&frequency<0.8]
#########DAPC#########################################################################################################################################################
#sortedcluster=ECOcleaned2A[,7:ncol(ECOcleaned2A)]
  sortedcluster=afterHWEtest
S1 <- nrow(sortedcluster)
EcoSNP <- vector(mode="list",S1)
for(i in 1:S1){
  Ecolistdata=as.integer(sortedcluster[i,])
  EcoSNP[[i]] <- new("SNPbin", Ecolistdata)
  print(i)
}
myPal=colorRampPalette(c("red","blue","black","green"))
Ecogen=new("genlight",EcoSNP,pop=rownames(sortedcluster))
grp=find.clusters(Ecogen,stat=c("BIC"))
png(file="findclusters.png",width=7,height=6,units="in",res=500)
plot(grp$Kstat, type="o", xlab="number of clusters (K)",
     ylab="BIC")
dev.off()
Ploidylist=DAPClist[,colnames(DAPClist)=='Ploidy']
Ploidylist[Ploidylist == 4]=1
Ploidylist[Ploidylist == 6]=2
Ploidylist[Ploidylist == 8]=3

grpsize=c(table(Ploidylist)[[1]],table(Ploidylist)[[2]],table(Ploidylist)[[3]])
Ploidylist=as.factor(Ploidylist)
attr(Ploidylist, "names")= as.character(c(1:165))
grp$grp=Ploidylist
grp$size=grpsize

dapcdata=dapc(Ecogen,grp$grp,n.pca=1,scale=F,var.contrib=T)
dapcdata$prior
dapcdata$pca.cent
dapcdata$assign
scatter(dapcdata, col=transp(myPal(6)), scree.da=FALSE,
        cell=1.5, cex=2, bg="white",cstar=0)
dapcdata$var.contr=100*dapcdata$var.contr
png(file="loadings.png",width=7,height=6,units="in",res=500)
loadingplot(dapcdata$var.contr, thres=0.15, ylab="Loadings(%)",cex.lab=0.8, xlab=NA_integer_, main=NA_integer_,frame.plot=FALSE)
dev.off()
box(col = "white")
axis(2, cex.axis = 1.5)



temp <- optim.a.score(dapcdata)

PPtable=cbind(dapcdata$posterior,dapcdata$assign,DAPClist)

write.csv(PPtable,file = "test7imputedDAPCPtable(PC=80,k=2,n.pca=1.new1).csv")

compoplot(dapcdata,col=c("red","blue","black","green"),lab=F, lwd=3,txt.leg=paste("group", 1:4), ncol=2)
colnames(PCAdata)[2804]
#########allele frequencies#########################################################################################################################################################
library("hierfstat")
ECOgenindUNIP=df2genind(converted, sep = "/", ind.names = newDAPClist[,2], ploidy =2)
strata(ECOgenindUNIP) <- my_strata
setPop(ECOgenindUNIP) =~Struclades
freq2962 <- tab(genind2genpop(ECOgenindUNIP[loc=c("MockRefGenome_11428639_G")]),freq=TRUE)
freq2963 <- tab(genind2genpop(ECOgenindUNIP[loc=c(colnames(PCAdata[2963]))]),freq=TRUE)
freq4972 <- tab(genind2genpop(ECOgenindUNIP[loc=c("MockRefGenome_40416897_G")]),freq=TRUE)
freq3143 <- tab(genind2genpop(ECOgenindUNIP[loc=c(colnames(PCAdata[3143]))]),freq=TRUE)
ECOgenpop=genind2genpop(ECOgenindUNIP)
colnames(PCAdata[4972])

genet.dist(ECOgenindUNIP, method = "WC84")
div <- summary(ECOgenindUNIP)
div
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
bartlett.test(list(div$Hexp, div$Hobs))


grp <- find.clusters(ECOgenindUNIP, max.n.clust = 100, n.pca = 100, choose.n.clust = T) 
names(grp)
dapc1 <- dapc(ECOgenindUNIP, grp$grp, n.pca = 30, n.da = 6) 
scatter(dapc1)

dev.off()
png("SNP#4972.png", width=5, height=7.2, units="in", res=500)
par(mfrow=c(1,1), mar=c(3.1,4.3,3.1,2.2),las=0)
matplot(freq4972, pch=NA_integer_,type="b",lwd=1.5,
        col=c("black","blue","red"),
        xlab="Ecotypes",ylab="allele frequency", xaxt="n",
        cex=1, main="SNP #4972")
text(freq4972[,1],label="cc",cex=1.2)
text(freq4972[,2],label="cg",cex=1.2)
text(freq4972[,3],label="gg",cex=1.2)
axis(side=1,at=1:3,lab=c("Central Mid-west and East","South Mid-west","North Mid-west"))
dev.off()

png("SNP#2962.png", width=5, height=7.2, units="in", res=500)
matplot(freq2962, pch=NA_integer_,type="b",lwd=1.5,
        col=c("black","blue","red"),
        xlab="Ecotypes",ylab="allele frequency", xaxt="n",
        cex=1, main="SNP #2962")
text(freq2962[,1],label="aa",cex=1.2)
text(freq2962[,2],label="ag",cex=1.2)
text(freq2962[,3],label="gg",cex=1.2)
axis(side=1,at=1:3,lab=c("Central Mid-west and East","South Mid-west","North Mid-west"))
dev.off()

png("SNP#2963.png", width=5, height=7.2, units="in", res=500)
matplot(freq2963, pch=NA_integer_,type="b",lwd=1.5,
        col=c("black","blue","red"),
        xlab="Ecotypes",ylab="allele frequency", xaxt="n",
        cex=1, main="SNP #2963")
text(freq2963[,1],label="aa",cex=1.2)
text(freq2963[,2],label="ac",cex=1.2)
text(freq2963[,3],label="cc",cex=1.2)
axis(side=1,at=1:3,lab=c("Central Mid-west and East","South Mid-west","North Mid-west"))
dev.off()

png("SNP#3143.png", width=5, height=7.2, units="in", res=500)
matplot(freq3143, pch=NA_integer_,type="b",lwd=1.5,
        col=c("black","blue","red"),
        xlab="Ecotypes",ylab="allele frequency", xaxt="n",
        cex=1, main="SNP #3143")
text(freq3143[,1],label="gg",cex=1.2)
text(freq3143[,2],label="ga",cex=1.2)
text(freq3143[,3],label="aa",cex=1.2)
axis(side=1,at=1:3,lab=c("Central Mid-west and East","South Mid-west","North Mid-west"))
dev.off()
#########PCA#########################################################################################################################################################
library(ape)
afterHWEtest=ECOcleaned2A[,7:ncol(ECOcleaned2A)]
rownames(DAPClist)=NULL
newDAPClist=DAPClist
rownames(newDAPClist)=NULL
PCAdata=afterHWEtest
rownames(PCAdata)=NULL
pca_matrix = scale(PCAdata) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})    #sets NAs to mean (0)
eig = prcomp(pca_matrix_noNA, center = T) #Performs PCA on pca_matrix_noNA
summary(eig) # print variance accounted for 
loadings(eig) # pc loadings 
plot(eig,type="lines") # scree plot 
screeplot(eig, type="lines")
eigenvalues = (eig$sdev)^2
scores = eig$x
plot(eig$x[,1],eig$x[,2]) 
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
pc2_var = round((eigenvalues[2]/sum(eigenvalues)) * 100, 1)
pc3_var = round((eigenvalues[3]/sum(eigenvalues)) * 100, 1)
plot(scores[,1],scores[,2]) # Plot PC1 versus PC2
PCV1=scores[,1]
PCV2=scores[,2]
PCV3=scores[,3]
myPCAdata=cbind(newDAPClist,PCV1,PCV2,PCV3)
rownames(myPCAdata)=NULL
#write.csv(myPCAdata,file = "SNPPCAdata.csv")
#myPCAdata=as.matrix(file.choose())
pcoa_matrix = scale(PCAdata) #adjusts mean and stdev of each column to 1 and 0 respectively
pcoa_matrix_noNA = apply(pcoa_matrix, 2, function (x) {ifelse(is.na(x), 0, x)}) 
PCAD=dist.gene(pcoa_matrix_noNA , method = "pairwise", pairwise.deletion = FALSE,
               variance = FALSE)
eig2=pcoa(PCAD, correction="none", rn=NULL)
PCOV1=eig2$vectors[,1]
PCOV2=eig2$vectors[,2]
PCOV3=eig2$vectors[,3]
biplot(pcoa(PCAD, correction="none", rn=NULL), Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL)
myPCOAdata=cbind(newDAPClist,PCOV1,PCOV2,PCOV3)

#########Scatterplots###############################################################################
fastclade=read.csv(file.choose(),header=F)
library(splitstackshape)
fastcladenew=cSplit(fastclade, 1:ncol(fastclade), sep=" ", type.convert=FALSE)
fastcladenew$V1_1=as.numeric(fastcladenew$V1_1)
fastcladenew$V1_2=as.numeric(fastcladenew$V1_2)
colnames(fastcladenew)=c("fast1","fast2","deme")
fastcladenew$deme[fastcladenew[,1]>0.5]=1
fastcladenew$deme[fastcladenew[,1]<0.5]=2

#myPCAdata$FastClade=as.factor(fastclade$FastClade)
#myPCAdata$Ploidy=as.factor(myPCAdata$Ploidy)

mapdata=cbind(myPCAdata,fastcladenew[,1:3])
mapdata[,8]=as.character(mapdata[,8])
mapdata[,9]=as.character(mapdata[,9])
mapdata$deme=as.factor(mapdata$deme)
#mapdata[50,8]='KS'
#mapdata[51,8]='KS'
#mapdata[52,8]='KS'
#mapdata[53,8]='KS'

mapdata[96,8]='Red River'
mapdata[96,9]='Red River'
gall <- subset(mapdata, State == "Red River")
gall_1 <- subset(mapdata, State2 == "NY2")
gall_2 <- subset(mapdata, State2 == "NY1")
gall_1[1,9]='KST'
gall_2[1,9]='STP'
gallnew=rbind(gall_1,gall_2,gall)
gallnew$State2=factor(gallnew$State2, levels=c('KST','STP','Red River'), ordered = TRUE)
mapdata=mapdata[-c(1,2,96),]
group.col=c("Blue","Red","orange","darkslategrey","lightseagreen","darkgoldenrod3","green3","cyan","burlywood4","darkred","orangered","purple","ivory4","plum4","black","springgreen2","hotpink4","lightskyblue4","darkorange","tan3","orchid4","sienna4")
png(file="TotalPop1.png",width=12,height=7,units="in",res=300)
q=ggplot(mapdata,aes(x=PCV1, y=PCV2))+
  geom_point(size=5,stroke = 0.5, aes(fill =deme,shape=as.factor(Ploidy)))+theme_bw()+xlab(paste("PC1(",pc1_var,"%)",sep=""))+ylab(paste("PC2(",pc2_var,"%)",sep=""))+
  geom_text(aes(label=State),hjust=0, vjust=-1.2)+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
    theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(1,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())+
  guides(shape = guide_legend(order=2))+
  geom_point(data=gallnew,aes(x=gallnew$PCV1,y=gallnew$PCV2,color=State2),fill=c("orange","black","black"),shape=c(18,7,9),size=c(8,5,6),stroke = 2.3) +
  scale_fill_manual(name = "Demes",values = c("Blue","Red"), labels=c("Deme1", "Deme2"))+
  guides(fill = guide_legend(override.aes = list(shape=18,color =c("Red","Blue"),size=8),order=1))+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(shape = guide_legend(override.aes = list(size=5,stroke=1.5),order=2))
  q=q+scale_color_manual(name="Released varieties",values =c("black","black","orange"))+
    guides(color = guide_legend(override.aes = list(shape=c(7,9,18),size=c(6,5,8),color =c("black","black","orange"))),order=3)
  q 
dev.off()
# scale_shape_manual(name = "F1&Parents",values = c(18,18,18,18,18,18,18,25,22,4,8,22,0,24,"#","%",9,3))
png(file="TotalPop2.png",width=12,height=7,units="in",res=300)
ggplot(mapdata, aes(x=PCV1, y=PCV2,fill =100*fast1,shape=as.factor(Ploidy)))+
  geom_point(size=5,stroke = 0.5)+theme_bw()+xlab(paste("PC1(",pc1_var,"%)",sep=""))+ylab(paste("PC2(",pc2_var,"%)",sep=""))+
  geom_text(aes(label=State),hjust=0, vjust=-1.2)+
    theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(0,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())+
  geom_point(data=gallnew,aes(x=gallnew$PCV1,y=gallnew$PCV2,color=State2),shape=c(7,9,18),size=c(5,6,8),stroke = 2.3)+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(fill = guide_colorbar(label.theme = element_text(angle=0,size=12,
                                                          face="bold"),order=1))+
  guides(shape = guide_legend(override.aes = list(size=5,stroke=1.5),order=2))+
  scale_fill_gradient(name = "%Posterior probability",low="red",high="blue",limits=c(0,100))+
  scale_colour_manual(name="Released varieties",values =c("black","black","orange"))+
  guides(colour = guide_legend(override.aes = list(shape=c(7,9,18),size=c(6,5,8),color =c("black","black","orange"))))
dev.off()


png(file="TotalPop2.png",width=12,height=7,units="in",res=300)
ggplot(mapdata, aes(x=PCV1, y=PCV2,fill =100*fast1,shape=as.factor(Ploidy)))+
  geom_point(size=5,stroke = 0.5)+theme_bw()+xlab(paste("PC1(",pc1_var,"%)",sep=""))+ylab(paste("PC2(",pc2_var,"%)",sep=""))+
  geom_text(aes(label=State),hjust=0, vjust=-1.2)+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(0,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())+
  geom_point(data=gallnew,aes(x=gallnew$PCV1,y=gallnew$PCV2,color=State2),shape=c(18,7,9),size=c(8,5,6),stroke = 2.3)+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(fill = guide_colorbar(label.theme = element_text(angle=0,size=12,
                                                          face="bold"),order=1))+
  guides(shape = guide_legend(override.aes = list(size=5,stroke=1.5),order=2))+
  scale_fill_gradient(name = "%Membership",low="red",high="blue",limits=c(0,100))+
  scale_colour_manual(name="Released varieties",values =c("black","black","orange"))+
  guides(colour = guide_legend(override.aes = list(shape=c(7,9,18),size=c(6,5,8),color =c("black","black","orange"))),order=3)
dev.off()

#3Dplots
newDAPClist=cbind(newDAPClist,PCV1,PCV2,PCV3)
my3ddataset=newDAPClist
my3ddataset[is.na(my3ddataset$mapcode),]$mapcode="Unknown"
my3ddataset[my3ddataset$mapcode=="Unknown",]$mapcode2="Unknown"


my3ddataset$groupcol=get_colors(my3ddataset$State, group.col)
open3d()
plot3d(my3ddataset[,17], my3ddataset[,18], my3ddataset[,19], type="n", box = FALSE, xlab="", ylab="", zlab="")
text3d(my3ddataset[,c(17,18,19)],text=my3ddataset$State,cex=1.8,fontweight='bold',family='sans')
rgl.postscript("3dplot.pdf", fmt="pdf")



x=my3ddataset[,16]
y=my3ddataset[,17]
z=my3ddataset[,18]


rgl_add_axes(my3ddataset[,16],my3ddataset[,17],my3ddataset[,18], show.bbox = TRUE)
#2Dlegend
png(file="legend.png",width=10,height=7,units="in",res=300)
plot.new()
legobject=ddply(my3ddataset, ~ mapcode, summarise, groupcol= min(groupcol))
legend("topright", legend =legobject$mapcode, pch = 15, pt.cex=2,col = legobject$groupcol, cex=1, inset=c(0.02))

dev.off()

M <- par3d("userMatrix")
if (!rgl.useNULL())
  play3d( par3dinterp(time = (0:2)*0.5, userMatrix = list(M,rotate3d(M, pi/0.2, 1, 1, 1),
                                                          rotate3d(M, pi/10, 0, 0, 0) ) ), 
          duration =10)
## Not run:
movie3d( spin3d(), duration = 5,dir = 'C:/Users/jiaguo3.UOFI/Box Sync/Box Sync/Jia/UNEAK codes/PNG')
## End(Not run)

# add legend
legend("topright", legend = my3ddataset$EcoregionIIINo, pch = 16, col = my3ddataset$EcoregionIIINo, cex=1)

# capture snapshot
snapshot3d(filename = '3dplot2.png', fmt = 'png')
library(plotly)
# volcano is a numeric matrix that ships with R
plot_ly(z = volcano, type = "surface")

#########AMOVA and fixationIndex#########################################################################################################################################################
library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("poppr")
library("PopGenReport")
library("hierfstat")
library(diveRsity)

newploidy=newDAPClist[,10]
newploidy[newploidy== 6] <- 8
myclades=read.csv(file.choose())

converted=PCAdata# x is a genlight object
converted[converted == 0] <- "1/1" # homozygote reference
converted[converted == 1] <- "1/2" # heterozygote
converted[converted == 2] <- "2/2" # homozygote alternate
ECOgenind <- df2genind(converted, sep = "/", ind.names = newDAPClist[,2], ploidy = newDAPClist[,10])
ECOgenind
my_strata <- data.frame(Population=newDAPClist[,9],State=newDAPClist[,8], Subsample=newDAPClist[,12],ploidy=newDAPClist[,10], fastclades=fastcladenew[,3])
strata(ECOgenind) <- my_strata
setPop(ECOgenind)=~Population
ECOgenind
ECOgenindnew=ECOgenind
ECOgenindnew$ploidy=as.integer(rep(2:2,len=96))
ECOgeninddeme1<- ECOgenindnew[ECOgenind$strata[,5] == "2", ]
ECOgeninddeme2<- ECOgenindnew[ECOgenind$strata[,5] == "1", ]

ECOgeninddeme1$ploidy=as.integer(rep(2:2,len=96))
ECOgeninddeme2$ploidy=as.integer(rep(2:2,len=96))

diff_stats(ECOgeninddeme1)
genetdistwc84=genet.dist(ECOgenindnew, method = "WC84")

setPop(ECOgenindnew) =~Population
deme1basic=basic.stats(ECOgeninddeme1)
deme2basic=basic.stats(ECOgeninddeme2)

deme1basic$overall
deme2basic$overall

toto <- summary(ECOgenindnew)
allel.rich(ECOgenindnew,min.alleles=NULL)
basic.stats(ECOgenindnew,digits=4)
fstat(ECOgenindnew)
data(genpopfile)
wc(ECOgenindnew)

str(toto)
mean(toto$Hobs) 
mean(toto$Hexp)

totodeme1= summary(ECOgeninddeme1)
setPop(ECOgeninddeme1) =~State
basic.stats(ECOgeninddeme1,digits=4)
genetdistwc84deme1=genet.dist(ECOgeninddeme1, method = "WC84")
wc(ECOgeninddeme1)
beep()
wc(ECOgeninddeme2)
beep()
mean(totodeme1$Hobs)
mean(totodeme1$Hexp)
write.csv(as.matrix(genetdistwc84deme1),file ="genetdistwc84deme1.csv")

totodeme2= summary(ECOgeninddeme2)
setPop(ECOgeninddeme2) =~State
basic.stats(ECOgeninddeme2,digits=4)
genetdistwc84deme2=genet.dist(ECOgeninddeme2, method = "WC84")
beep()
mean(totodeme2$Hobs)
mean(totodeme2$Hexp)
write.csv(as.matrix(genetdistwc84deme2),file ="genetdistwc84deme2.csv")


allFst=pairwise.fst(ECOgenindnew)
str(allFst)
write.csv(as.matrix(allFst),file ="allFst.csv")
write.csv(as.matrix(genetdistwc84),file ="genetdistwc84.csv")


allGfst=pairwise_Gst_Nei(ECOgenindnew, linearized = T)
deme1fst=pairwise_Gst_Nei(ECOgeninddeme1, linearized = T)
deme2fst=pairwise_Gst_Nei(ECOgeninddeme2, linearized = T)
deme1fstHedrick=pairwise_Gst_Hedrick(ECOgeninddeme1, linearized = FALSE)
pairwise_D(ECOgenind, linearized = FALSE, hsht_mean = "harmonic")
bs <- chao_bootstrap(ECOgenind, nreps = 100)
summarise_bootstrap(bs, Gst_Nei)     # for Nei's Gst
summarise_bootstrap(bs, Gst_Hedrick) # for Hedrick's Gst
summarise_bootstrap(bs, D_Jost)      # for Jost's D

testECOgenid=ECOgenind
amova.result=poppr.amova(testECOgenid, ~ploidy/Population, filter = TRUE, threshold = 0.1, clonecorrect = TRUE)
amova.result
Ecosignif=randtest(amova.result , nrepet = 999)
amova.result2<- poppr.amova(testECOgenid, ~Struclades/Population, filter = TRUE, method='pegas', nperm=10000,threshold = 0.1,clonecorrect=T)

library(radiator)
library(assigner)
tidygenind=tidy_genind(ECOgenindnew)
iteratedfst.pairwise <- fst_WC84(
  tidygenind, 
  pairwise = TRUE,
  ci = TRUE, 
  iteration.ci = 10000, 
  quantiles.ci = c(0.025,0.975),
  parallel.core = 8,
  verbose = TRUE
)
wc84iterate=iteratedfst.pairwise$pairwise.fst.ci.matrix
iteratedfst.pairwise$fst.plot
write.csv(as.matrix(wc84iterate),file ="wc84iterate.csv")

#########tree#########################################################################################################################################################

library(ape)
newDAPClist[,8]=as.character(newDAPClist[,8])
newDAPClist[50,8]='KS'
newDAPClist[51,8]='KS'
newDAPClist[52,8]='KS'
newDAPClist[53,8]='KS'
treedata=afterHWEtest
rownames(treedata)=paste(newDAPClist[,9],"-",newDAPClist[,12])
treefit<- hclust(dist(as.matrix(treedata)), members = NULL)

nclus=3
color=c('red','blue','black',"purple","orange")
color_list=rep(color,nclus/length(color))
clus=cutree(treefit,nclus)


plot(as.phylo(treefit),type='radial',tip.color=color_list[clus],label.offset=0.2,no.margin=TRUE, cex=0.70, show.node.label = TRUE)


plot.fan <- function(hc, nclus=9) {
  palette <- c("blue",'red','blue','red',"purple","gold","purple","purple","red")[1:nclus]
  clus    <-cutree(hc,nclus)
  X <- as.phylo(hc)
  edge.clus <- sapply(1:nclus,function(i)max(which(X$edge[,2] %in% which(clus==i))))
  order     <- order(edge.clus)
  edge.clus <- c(min(edge.clus),diff(sort(edge.clus)))
  edge.clus <- rep(order,edge.clus)
  plot(X,type = "radial",
       tip.color=palette[clus],edge.color=palette[edge.clus],
       show.tip.label = T,no.margin=TRUE, cex=0.70)
}
png("Tree.png", width=10, height=8, units="in", res=1000)
plot.fan(treefit ,2)
dev.off()
plot.phylo2(treefit,show.tip.label = T,use.edge.length = T, no.margin = T, cex = 0.55)

cladedata=as.matrix(clus)
newDAPClist$hclus=as.factor(cladedata[,1])

#########genetic distance#########################################################################################################################################################

library(cluster)
library(reshape)
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}
GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

#ifnotimputed: missingno(testECOgenid, type = "mean")
#ECO_dist=diss.dist(testECOgenid)
#amova(Popdist~ploidy, data=strata(ECOgenind),nperm=100)
geodistancetemp=cbind(myPCAdata,fastcladenew[,1:3])
Geneticdistancedata=PCAdata
newDAPClist$FastClade=as.factor(fastclade$FastClade)
newDAPClist$DPACClade=as.factor(dapcdata$assign)
rownames(Geneticdistancedata)=paste(geodistancetemp[,18],"-",geodistancetemp[,9],"-",geodistancetemp[,10])
#rownames(Geneticdistancedata)=paste(my_strata[,4],"-",my_strata[,5],"-",my_strata[,1],my_strata[,2])
Popdist=as.dist(daisy(Geneticdistancedata,metric = c("gower")))
Popdistlist=melt(as.matrix(Popdist))[melt(upper.tri(as.matrix(Popdist)))$value,]
#correlation

geodistance=cbind(paste(geodistancetemp[,18],"-",geodistancetemp[,9],"-",geodistancetemp[,10]),geodistancetemp[,c(6,7)])
colnames(geodistance)=c("index","lat","lon")
geodistanceM=as.dist(round(GeoDistanceInMetresMatrix(geodistance) / 1000))
geodistlist=melt(as.matrix(geodistanceM))[melt(upper.tri(as.matrix(geodistanceM)))$value,]
#mantel.rtest(geodistanceM,Popdist)
#plot(r1 <- mantel.rtest(geodistanceM,Popdist), main = "Mantel's test")
#r1
distcorr=cbind(Popdistlist,geodistlist)
distcorrplot=distcorr[,-c(4,5)]
distcorrplot$X1=as.character(distcorrplot$X1)
distcorrplot$X2=as.character(distcorrplot$X2)
distcorrnew=cbind(do.call(rbind, strsplit(distcorrplot[,1], '\\-')),do.call(rbind, strsplit(distcorrplot[,2], '\\-')),distcorrplot[,3:4])
colnames(distcorrnew)=c("Geno1_deme","Geno1_Origin","Geno1_Ploidy","Geno2_deme","Geno2_Origin","Geno2_Ploidy","Genetic_distance","Geographicl_distance")
distcorrnew1=distcorrnew[distcorrnew[,1]==distcorrnew[,4],]
tetradistcorr=distcorrnew1[distcorrnew1[,3]==distcorrnew1[,6],]
plot(distcorrnew1[,8],distcorrnew1[,7], xlab="geograhical distance",ylab="genetic distance",main="Scatterplot Example", pch=19)
ggplot(tetradistcorr, aes(x=Geographicl_distance, y=Genetic_distance, group=Geno1_deme)) +geom_point(aes(shape=Geno1_Ploidy, color=Geno1_deme))

write.csv(distcorr,file = "distcorr.csv")
boxplot()
#########heterozygosity#####
chtest=as.data.frame(as.matrix(cbind(mylist[,10],PCAdata)))
chtest=chtest[-c(11,12,13,22),]
fishertable= c() 
for (i in 2:ncol(PCAdata)) {
  print(i)
  
  Fourhomo_ref_sum = sum(chtest[chtest[,1]==4,i] == 0, na.rm = T)
  Fourhetero_sum = sum(chtest[chtest[,1]==4,i] == 1, na.rm = T)
  Fourhomo_alt_sum = sum(chtest[chtest[,1]==4,i] == 2,na.rm=T)

  Eighthomo_ref_sum = sum(chtest[chtest[,1]==8,i] == 0, na.rm = T)
  Eighthetero_sum = sum(chtest[chtest[,1]==8,i] == 1, na.rm = T)
  Eighthomo_alt_sum = sum(chtest[chtest[,1]==8,i] == 2,na.rm=T)
  
  Eightobserved_geno_freqs = t(as.matrix(c(Eighthomo_ref_sum, Eighthetero_sum, Eighthomo_alt_sum)))
  Fourobserved_geno_freqs = t(as.matrix(c(Fourhomo_ref_sum, Fourhetero_sum, Fourhomo_alt_sum)))
  
  homo_ref_sum = sum(chtest[,i] == 0, na.rm = T)
  hetero_sum = sum(chtest[,i] == 1, na.rm = T)
  homo_alt_sum = sum(chtest[,i] == 2,na.rm=T)

  sample_size = sum(homo_ref_sum, hetero_sum, homo_alt_sum)

  ref_freq = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
  alt_freq = ((homo_alt_sum * 2) + hetero_sum)/(sample_size * 2)

  expected_hetero = 2*ref_freq*alt_freq
  observed_hetero = hetero_sum/sample_size
  
  freqtable=as.matrix(rbind(Eightobserved_geno_freqs,Fourobserved_geno_freqs))
  fishertable[,i-1]=fisher.test(freqtable)[1]
}
table(fishertable<0.05)
library(plyr)
library(multtest)
fishertablenew=as.matrix(t(as.data.frame(fishertable)))
rownames(fishertablenew)=NULL
fishertablenew=cbind(fishertablenew,p.adjust(fishertablenew[,1], method='BH'))
table(fishertablenew[,2]<0.01)

res<- mt.rawp2adjp(fishertablenew[,1], proc=c("BH"),alpha=0.001,na.rm = FALSE)
adjp <- res$adjp[order(res$index), ]
dim(na.omit(adjp))
table(na.omit(adjp[,2]<0.001))
pvalues[,4]=adjp[,2]

heterdata=as.data.frame(as.matrix(cbind(fastcladenew[,3],PCAdata)))
deme1pop=heterdata[heterdata[,1]==2,]
deme2pop=heterdata[heterdata[,1]==1,]
deme1obhetero=c()
deme1exphetero=c()
deme1fsh=c()
pb <- txtProgressBar(min = 5, max = ncol(chtest), style = 3)
for (i in 2:ncol(chtest)) {
  homo_ref_sum = sum(deme1pop[,i] == 0, na.rm = T)
  hetero_sum = sum(deme1pop[,i] == 1, na.rm = T)
  homo_alt_sum = sum(deme1pop[,i] == 2,na.rm=T)

  sample_size = sum(homo_ref_sum, hetero_sum, homo_alt_sum)
  ref_freq = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
  alt_freq = ((homo_alt_sum * 2) + hetero_sum)/(sample_size * 2)
  
  expected_homo_ref = ref_freq^2
  expected_hetero = 2*ref_freq*alt_freq
  expected_homo_alt = alt_freq^2
  
  deme1obhetero[i]= hetero_sum/sample_size
  deme1exphetero[i]=expected_hetero
  deme1fsh[i]=(deme1exphetero[i]-deme1obhetero[i])/deme1exphetero[i]
  
  setTxtProgressBar(pb, i+5)
}

boxplot(deme1obhetero)
sd(deme1obhetero)
mean(deme1obhetero)
boxplot(deme1exphetero)
boxplot(deme1fsh)

deme2obhetero=c()
deme2exphetero=c()
deme2fsh=c()
pb <- txtProgressBar(min = 5, max = ncol(chtest), style = 3)
for (i in 2:ncol(chtest)) {
  homo_ref_sum = sum(deme2pop[,i] == 0, na.rm = T)
  hetero_sum = sum(deme2pop[,i] == 1, na.rm = T)
  homo_alt_sum = sum(deme2pop[,i] == 2,na.rm=T)
  
  sample_size = sum(homo_ref_sum, hetero_sum, homo_alt_sum)
  ref_freq = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
  alt_freq = ((homo_alt_sum * 2) + hetero_sum)/(sample_size * 2)
  
  expected_homo_ref = ref_freq^2
  expected_hetero = 2*ref_freq*alt_freq
  expected_homo_alt = alt_freq^2
  
  deme2obhetero[i]= hetero_sum/sample_size
  deme2exphetero[i]=expected_hetero
  deme2fsh[i]=(deme1exphetero[i]-deme1obhetero[i])/deme1exphetero[i]
  
  setTxtProgressBar(pb, i+5)
}

boxplot(deme2obhetero)
boxplot(deme2exphetero)
boxplot(deme2fsh)