library(adegenet)
library(scrime)
library(ggplot2)
library(rgl)
library(animation)
library(ape)
library(plyr)
library(reshape2)
library(poppr)
library(ade4)
hapMap2genlight <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)]
  samples <- scan(file, what = character(), nlines = 1)[-(1:11)]
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
##########################################################################################################################################################
#input data
mydata <- hapMap2genlight(file.choose())
mydata2=as.matrix(mydata)
rownames(mydata2)
dim(mydata2)

ECOdata=mydata2[c(1:4,6:212),]
ECOdata=mydata2[c(77:285),]

dim(ECOdata)
datasetname=as.matrix(rownames(ECOdata))

mylist=read.csv(file.choose())
rownames(ECOdata)=mylist[,1]


str(mydata2)
NaCountingforrow=matrix(NA,nrow(mydata2),3)

for (i in 1:nrow(mydata2)){
  NaCountingforrow[i,1]=sum(is.na(mydata2[i,]))
  NaCountingforrow[i,2]=rownames(mydata2)[i]
}
write.csv(NaCountingforrow,file = "countingdata.csv")
str(NaCountingforrow)
table(mydata2[255,]==mydata2[254,])
##########################################################################################################################################################
#quality control
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
for (i in 1:ncol(newECOdata)) {
  print(i)
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
  }

sum(as.numeric(matchingdata[134,-1]==matchingdata[154,-1]),na.rm = T)/length(which(!is.na(matchingdata[154,-1])))
#Allele frequency test
frequency=c()
for (i in 1:ncol(mydata2)){
  frequency[i]=round(sum(mydata2[,i])/90,2)  
  print(i)
}
newmydata2=mydata2[,frequency>0.2&frequency<0.8]
##########################################################################################################################################################
#ouputSNPdata
output=afterHWEtest
output[is.na(output)]=-9
#my3ddataset$NOEcoregion=as.integer(my3ddataset[,colnames(my3ddataset)=='mapcode2'])
output2=cbind(as.matrix(my3ddataset$ID),as.matrix(my3ddataset$mapcode2),output)
write.table(output2, file="newECO2.txt",sep = "\t")
##########################################################################################################################################################
#AMOVA
Vardata=cbind(mylist$PopNo)






##########################################################################################################################################################
#PCA
pca_matrix = scale(afterHWEtest) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})    #sets NAs to mean (0)
eig = prcomp(pca_matrix_noNA, center = T) #Performs PCA on pca_matrix_noNA
summary(eig) # print variance accounted for 
loadings(eig) # pc loadings 
plot(eig,type="lines") # scree plot 
screeplot(eig, type="lines")
eigenvalues = (eig$sdev)^2
scores = eig$x 
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
pc2_var = round((eigenvalues[2]/sum(eigenvalues)) * 100, 1)
pc3_var = round((eigenvalues[3]/sum(eigenvalues)) * 100, 1)
plot(scores[,1],scores[,2]) # Plot PC1 versus PC2
PCV1=scores[,1]
PCV2=scores[,2]
PCV3=scores[,3]
myPCAdata=cbind(newmylist,PCV1,PCV2,PCV3)
rownames(myPCAdata)=NULL
#write.csv(myPCAdata,file = "SNPPCAdata.csv")
#myPCAdata=as.matrix(file.choose())
#Scatterplots
################################################################################
#2Dplots
group.col=c("yellowgreen","magenta4","blue","darkslategrey","lightseagreen","darkgoldenrod3","green3","cyan","burlywood4","darkred","orangered","purple","ivory4","plum4","springgreen2","hotpink4","black","lightskyblue4","darkorange","tan3","orchid4","sienna4")
png(file="ECOdata.png",width=10,height=7,units="in",res=300)
ggplot(myPCAdata, aes(x=PCV1, y=PCV2,colour =EcoregionIII,group=EcoregionIII))+
  geom_point(size=3)+theme_bw()+xlab(paste("PC1(",pc1_var,"%)",sep=""))+ylab(paste("PC2(",pc2_var,"%)",sep=""))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"))+
  scale_colour_manual(name = "EcoregionIII",values = group.col)
dev.off()
# scale_shape_manual(name = "F1&Parents",values = c(18,18,18,18,18,18,18,25,22,4,8,22,0,24,"#","%",9,3))

#3Dplots
my3ddataset=myPCAdata
my3ddataset[is.na(my3ddataset$mapcode),]$mapcode="Unknown"
my3ddataset[my3ddataset$mapcode=="Unknown",]$mapcode2="Unknown"


my3ddataset$groupcol=get_colors(my3ddataset$mapcode, group.col)
 open3d()
 text3d(my3ddataset[,c(12,13,14)],col=my3ddataset$groupcol,text=my3ddataset$State,cex=1.2,type='s',fontweight='bold')
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
#########################################################################################################################################################
#imputation
newmydata2<- knncatimputeLarge(t(mydata2))
#second allele
myFreq <- glMean(mydata)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)
#first allele
myFreq <- glMean(mydata)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
#missing allele
temp <- density(glNA(mydata), bw=10)
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
     xlim=c(0,1701))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(glNA(flu), rep(0, nLoc(flu)), pch="|", col="blue")
#missing allele for each family
missings=NA.posi(mydata)
##########################################################################################################################################################
#tree
treedata=afterHWEtest
rownames(treedata)=newmylist[rownames(treedata),8]
treefit<- hclust(dist(as.matrix(treedata)), method = "complete", members = NULL)

nclus=4
color=c('red','blue','black',"purple","orange")
color_list=rep(color,nclus/length(color))
clus=cutree(treefit,nclus)


plot(as.phylo(treefit),type='radial',tip.color=color_list[clus],label.offset=0.2,no.margin=TRUE, cex=0.70, show.node.label = TRUE)


plot.fan <- function(hc, nclus=4) {
  palette <- c('red','blue','black','orange',"purple","cyan","hotpink4","darkslategrey")[1:nclus]
  clus    <-cutree(hc,nclus)
  X <- as.phylo(hc)
  edge.clus <- sapply(1:nclus,function(i)max(which(X$edge[,2] %in% which(clus==i))))
  order     <- order(edge.clus)
  edge.clus <- c(min(edge.clus),diff(sort(edge.clus)))
  edge.clus <- rep(order,edge.clus)
  plot(X,type='radial',
       tip.color=palette[clus],edge.color=palette[edge.clus],
       label.offset=0.02,no.margin=TRUE, cex=0.7)  
}

plot.fan(treefit ,4)


cladedata=as.matrix(clus)
my3ddataset$clade=as.factor(cladedata[,1])

##########################################################################################################################################################
#DAPC
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

dapcdata=dapc(Ecogen,grp$grp,n.pca=1,scale=FALSE)
scatter(dapcdata, col=transp(myPal(6)), scree.da=FALSE,
        cell=1.5, cex=2, bg="white",cstar=0)
loadingplot(dapcdata$var.contr, thres=1.3e-3)
temp <- optim.a.score(dapcdata)

PPtable=cbind(dapcdata$posterior,dapcdata$assign,newmylist)

write.csv(PPtable,file = "CombindPPtable(PC=100,k=4,n.pca=1.new).csv")

compoplot(dapcdata,col=c("red","blue","black","green"),lab=F, lwd=3,txt.leg=paste("group", 1:4), ncol=2)

##########################################################################################################################################################
#geneticVSphyical distance
dist.genpop(mydata, method = 1, diag = FALSE, upper = FALSE)
