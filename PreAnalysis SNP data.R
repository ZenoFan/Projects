#####################################convert to genlight##############################
hapMap2genlight2 <- function(file){
  require(adegenet)
  hapmapraw1 <- read.csv(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE,check.names = F)[,-(2:10)]
  #monolist=c("A","C","G","T")
  delelist=c("A/-/+","A/+/-","-/+/A","-/A/+","+/-/A",
             "C/-/+","C/+/-","-/+/C","-/C/+","+/-/C",
             "G/-/+","G/+/-","-/+/G","-/G/+","+/-/G",
             "T/-/+","T/+/-","-/+/T","-/T/+","+/-/T",
             "-/+","+/-")
  #hapmapraw2 <- hapmapraw1[!hapmapraw1$alleles%in%monolist, ]
  #table(hapmapraw1$alleles%in%monolist)
  hapmap <- rbind(hapmapraw1[!hapmapraw1$alleles%in%delelist, ],hapmapraw1[hapmapraw1$alleles%in%delelist, ])
  delecount=table(hapmapraw1$alleles%in%delelist)["TRUE"]
  
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
  #adele=s
  #cdele=s
  #gdele=s
  #tdele=s
  adeleN=s
  cdeleN=s
  gdeleN=s
  tdeleN=s
  dele=s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  #names(adele) <- c("-","0","A","N")
  #names(cdele) <- c("-","0","C","N")
  #names(gdele) <- c("-","0","G","N")
  #names(tdele) <- c("-","0","T","N")
  names(dele) <- c("-","0","Dummy","N")
  
  names(adeleN) <- c("-","0","A","N")
  names(cdeleN) <- c("-","0","C","N")
  names(gdeleN) <- c("-","0","G","N")
  names(tdeleN) <- c("-","0","T","N")
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt,
               #adele,adele,cdele,cdele,gdele,gdele,tdele,tdele,dele,dele,
               adeleN,adeleN,adeleN,adeleN,adeleN,
               cdeleN,cdeleN,cdeleN,cdeleN,cdeleN,
               gdeleN,gdeleN,gdeleN,gdeleN,gdeleN,
               tdeleN,tdeleN,tdeleN,tdeleN,tdeleN,
               dele,dele)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G",
                   "A/-/+","A/+/-","-/+/A","-/A/+","+/-/A",
                   "C/-/+","C/+/-","-/+/C","-/C/+","+/-/C",
                   "G/-/+","G/+/-","-/+/G","-/G/+","+/-/G",
                   "T/-/+","T/+/-","-/+/T","-/T/+","+/-/T",
                   "-/+","+/-")
  
  #"A/-","-/A","C/-","-/C","G/-","-/G","T/-","-/T","-/+","+/-"
  # Pull out and convert genotypes
  S <- length(samples)
  SBlist <- vector(mode="list",S)   # set up list of SNPbin objects
  
  total <- S
  # create progress bar
  pb <- txtProgressBar(min = 1, max = S, style = 3)

  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=T, USE.NAMES=FALSE)
    mygen[sapply(mygen , is.null)] <- NA
    # create SNPbin object for this individual
    SBlist[[i]] <- new("SNPbin", unlist(mygen))
    setTxtProgressBar(pb, i)
  }
  
  # make genlight object
  x <- new("genlight", SBlist,n.loc=length(loci),ind.names=samples,loc.names=loci)
output=list(x,delecount)
  return(output)
}
mygeno=hapMap2genlight2(file.choose())
###############################HWE filter############################################
mygenodata=as.data.frame(mygeno[1])
tempgenodata=mygenodata[,c((ncol(mygenodata)-as.numeric(mygeno[2])+1):ncol(mygenodata))]
HWEdata=mygenodata[,1:(ncol(mygenodata)-as.numeric(mygeno[2]))]
chisq_pval=c()
pb <- txtProgressBar(min = 1, max = ncol(HWEdata), style = 3)
for (i in 1:ncol(HWEdata)) {
  #print(i)
  #Calculate the genotype frequencies for each of the three possible genotypes.  Note that "na.rm = T" must be used to exclude NAs
  homo_ref_sum = sum(HWEdata[,i] == 0, na.rm = T)
  hetero_sum = sum(HWEdata[,i] == 1, na.rm = T)
  homo_alt_sum = sum(HWEdata[,i] == 2,na.rm=T)
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
  setTxtProgressBar(pb, i)
}
afterHWEgenodata=cbind(HWEdata[,chisq_pval>=0.05],tempgenodata)
###################################PCA######################################################
library(ggplot2)
pca_matrix = scale(mygenodata)#use mygenodata if HWE is not a issue
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)}) 
eig = prcomp(pca_matrix_noNA, center = T)
summary(eig)
loadings(eig)
plot(eig,type="lines")
eigenvalues = (eig$sdev)^2
scores = eig$x
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1)
pc2_var = round((eigenvalues[2]/sum(eigenvalues)) * 100, 1)
pc3_var = round((eigenvalues[3]/sum(eigenvalues)) * 100, 1)
plot(scores[,1],scores[,2])
PCV1=scores[,1]
PCV2=scores[,2]
PCV3=scores[,3]

PCAdata=as.data.frame(cbind(as.numeric(PCV1),as.numeric(PCV2),as.numeric(PCV3)))
PCAdata$Pop=as.vector(rownames(mygenodata))#use mygenodata if HWE is not a issue
png(file="PCAsumit.png",width=15,height=11,units="in",res=700)#output picture
ggplot(PCAdata,aes(x=PCV1, y=PCV2))+
  geom_point(size=5,stroke = 0.5)+theme_bw()+xlab(paste("PC1(",pc1_var,"%)",sep=""))+ylab(paste("PC2(",pc2_var,"%)",sep=""))+
  geom_text(aes(label=Pop),hjust=0, vjust=-1.2)+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(1,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())
dev.off()
################################kinship####################
library("scales")
#with EMMA method(load EMMA first)
foremmakinship=afterHWEgenodata#use mygenodata if HWE is not a issue
foremmakinship=rescale(as.matrix(foremmakinship), to=c(0,1))
k_matrix=emma.kinship(t(foremmakinship),"additive","all") 
rownames(k_matrix)=rownames(foremmakinship)
heatmap(as.matrix(K_matrix), scale='none')
#use kinship from TASSEL
tasselKm<- read.csv(file.choose(), header=F, row.names=1, sep="\t",
                       stringsAsFactors=FALSE,check.names=F)
heatmap(as.matrix(tasselKm), scale='none')
##############################LDdecaypre######################################################
hapmapinfo<- read.csv(file.choose(), header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,1:3]
hapmapinfo=hapmapinfo[rownames(hapmapinfo)%in%colnames(afterHWEgenodata),]
LDmap=hapmapinfo
LDmap[,1]=rownames(hapmapinfo)
colnames(LDmap)=c("Locus","LG","Position")
rownames(LDmap)=NULL
library(LDcorSV)
###LDdecay1A#########################################################
LDmap1A=LDmap[LDmap[,2]=="1A",]
LDdata1A=afterHWEgenodata[,colnames(afterHWEgenodata)%in%LDmap1A[,1]]
LD1A=LD.Measures(LDdata1A, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
#LD1A$distance=NA
#for (i in 1:nrow(LD1A)){
#  print(i)
 #  LD1A[i,4]=LDmap[LDmap[1]==as.matrix(LD1A)[i,1],3]
#}
LD1A$distance=NA
LD1A$distance=abs(as.numeric(gsub("S1A_", "", LD1A[,1]))-as.numeric(gsub("S1A_", "", LD1A[,2])))
png(file="LD1A.png",width=12.57703333,height=8,units="in",res=700)
plot(LD1A$distance,LD1A$r2)
dev.off()
lo <- loess(LD1A$r2~LD1A$distance)
plot(LD1A$distance,LD1A$r2)
lines(predict(lo), col='red', lwd=2)
###LDdecay2B#########################################################
LDmap2B=LDmap[LDmap[,2]=="2B",]
LDdata2B=afterHWEgenodata[,colnames(afterHWEgenodata)%in%LDmap2B[,1]]
LD2B=LD.Measures(LDdata2B, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
LD2B$distance=NA
LD2B$distance=abs(as.numeric(gsub("S2B_", "", LD2B[,1]))-as.numeric(gsub("S2B_", "", LD2B[,2])))

png(file="LD2B.png",width=12.57703333,height=8,units="in",res=700)
plot(LD2B$distance,LD2B$r2)
dev.off()

lo <- loess(LD2B$r2~LD2B$distance)
plot(LD2B$distance,LD2B$r2)
lines(predict(lo), col='red', lwd=2)
###LDdecay3D#########################################################
LDmap3D=LDmap[LDmap[,2]=="3D",]
LDdata3D=afterHWEgenodata[,colnames(afterHWEgenodata)%in%LDmap3D[,1]]
LD3D=LD.Measures(LDdata3D, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
LD3D$distance=NA
LD3D$distance=abs(as.numeric(gsub("S3D_", "", LD3D[,1]))-as.numeric(gsub("S3D_", "", LD3D[,2])))
png(file="LD3D.png",width=12.57703333,height=8,units="in",res=700)
plot(LD3D$distance,LD3D$r2)
dev.off()

lo <- loess(LD3D$r2~LD3D$distance)
plot(LD3D$distance,LD3D$r2)
lines(predict(lo), col='red', lwd=2)
###LDdecayA#########################################################
LDmapA=LDmap[substr(LDmap[,2],2,2)=="A",]
LDdataA=afterHWEgenodata[,colnames(afterHWEgenodata)%in%LDmapA[,1]]
LDA=LD.Measures(LDdataA, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
LDA$distance=NA
LDA$distance=abs(as.numeric(substr(LDA[,1],5,100))-as.numeric(substr(LDA[,2],5,100)))
png(file="LDA.png",width=12.57703333,height=8,units="in",res=700)
plot(LDA$distance,LDA$r2)
dev.off()

lo <- loess(LDA$r2~LDA$distance)
plot(LDA$distance,LDA$r2)
lines(predict(lo), col='red', lwd=2)
###LDdecayD#########################################################
LDmapD=LDmap[substr(LDmap[,2],2,2)=="D",]
LDdataD=afterHWEgenodata[,colnames(afterHWEgenodata)%in%LDmapD[,1]]
LDD=LD.Measures(LDdataD, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
LDD$distance=NA
LDD$distance=abs(as.numeric(substr(LDD[,1],5,100))-as.numeric(substr(LDD[,2],5,100)))
png(file="LDD.png",width=12.57703333,height=8,units="in",res=700)
plot(LDD$distance,LDD$r2)
dev.off()
###LDdecayB#########################################################
LDmapB=LDmap[substr(LDmap[,2],2,2)=="B",]
LDdataB=afterHWEgenodata[,colnames(afterHWEgenodata)%in%LDmapB[,1]]
LDB=LD.Measures(LDdataB, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
LDB$distance=NA
LDB$distance=abs(as.numeric(substr(LDB[,1],5,100))-as.numeric(substr(LDB[,2],5,100)))
png(file="LDB.png",width=12.57703333,height=8,units="in",res=700)
plot(LDB$distance,LDB$r2)
dev.off()
