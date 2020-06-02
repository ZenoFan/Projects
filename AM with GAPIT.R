# GAPIT - Genomic Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka, Feng Tian and You Tang
# Last update: September 15, 2015

#Install packages (Do this section only for new installation of R)
#-------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
library(BiocManager)
BiocManager::install()
#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#source("/Users/Zhiwu/Dropbox/Current/revolutionr/gapit/gapit_functions.txt")
#############################################################################################

is.even <- function(x) x %% 2 == 0
DipY  <- read.csv("C:/Users/Dipendra/Google Drive/Ph.D life/Analysis/Dipendra Analysis/Together_data/Citra_2018/Citra2018_lsmeans.csv", header=T, sep=",",
                  stringsAsFactors=FALSE,check.names = F)
sumitYlist  <- read.csv("D:/Dropbox (UFL)/wheat/Sumit_Jia/Previous/sumit/entry2genotype.csv", header=T, sep=",",
                    stringsAsFactors=FALSE,check.names = F)
DipG <- read.csv("C:/Users/Dipendra/Google Drive/Ph.D life/Analysis/Dipendra Analysis/Together_data/4_23_filtered.hmp.txt",header=F, sep="\t",
         stringsAsFactors=FALSE)

str(DipY)

DipY[,1]=gsub("Entry","",DipY[,1])
DipY[,1]=as.numeric(DipY[,1])

colnames(DipY)[1]="Taxa"

DipY[,1]=sumitYlist$GENOTYPE[match(DipY[,1], sumitYlist[,1])]

DipG[1,1:11]=c("rs","alleles","chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panel", "QCcode")
DipG=DipG[,DipG[1,-(1:11)]%in%DipY[,1]]

DipY=DipY[DipY[,1]%in%DipG[1,-(1:11)],]
target <- DipG[1,-(1:11)]
DipY=DipY[match(target, DipY$Taxa),]
#DipY=DipY[,-1]
rownames(DipY)=NULL

pca_matrix = scale(DipY[,c(4,5,6,20)]) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})    #sets NAs to mean (0)
eig = prcomp(pca_matrix_noNA, center = T) #Performs PCA on pca_matrix_noNA
eigenvalues = (eig$sdev)^2
scores = eig$x
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
plot(scores[,1],scores[,2]) # Plot PC1 versus PC2
PCV1=scores[,1]
DipYnew=cbind(DipY,PCV1)
colnames(DipYnew)[33]="BiomassPCA1"

pca_matrix = scale(DipY[,c(7,8,9)]) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})    #sets NAs to mean (0)
eig = prcomp(pca_matrix_noNA, center = T) #Performs PCA on pca_matrix_noNA
eigenvalues = (eig$sdev)^2
scores = eig$x
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
plot(scores[,1],scores[,2]) # Plot PC1 versus PC2
PCV1=scores[,1]
DipYnew=cbind(DipYnew,PCV1)
colnames(DipYnew)[34]="Reproductive partitioningPCA1"



#Using ECMLM by Li and et. al. (BMC Biology, 2014)
myGAPIT <- GAPIT(
  Y=DipY,
  G=DipG,
  kinship.cluster=c("average", "complete", "ward"),
  kinship.group=c("Mean", "Max"),
  group.from=200,
  group.to=1000000,
  group.by=10
)


Dipsnplist <- read.csv("D:/Dropbox (UFL)/wheat/Jia/2017 KASP study/report/4_23 snp list.csv", header=T, sep=",",
                 stringsAsFactors=FALSE)
DipG1=DipG[-1,]

Dipsnplist$allele=DipG1[match(Dipsnplist$SNP, DipG1$V1),2]
Dipsnplist$chr=DipG1[match(Dipsnplist$SNP, DipG1$V1),3]
Dipsnplist$seq=DipG1[match(Dipsnplist$SNP, DipG1$V1),12:253]
colnames(Dipsnplist$seq)=DipG[1,12:253]

write.csv(Dipsnplist ,"D:/Dropbox (UFL)/wheat/Jia/2017 KASP study/report/snp list with genotypes_4_23.csv")

selected_snplist <- read.csv("D:/Dropbox (UFL)/wheat/Jia/2017 KASP study/report/snp list with genotypes_4_23_selected.csv", header=T, sep=",",
                       stringsAsFactors=FALSE)
matchlist <- read.csv("D:/Dropbox (UFL)/wheat/Jia/2017 KASP study/report/match list.csv", header=T, sep=",",
                             stringsAsFactors=FALSE)

selected_snplist$picked=ifelse((selected_snplist$SNP)%in%(matchlist$snp),"picked","")
write.csv(selected_snplist ,"D:/Dropbox (UFL)/wheat/Jia/2017 KASP study/report/snp list with genotypes_4_23_selected_2.csv")
