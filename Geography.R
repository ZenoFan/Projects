library(sp)
library(rgdal) 
library(PBSmapping)
library(maptools)
library(raster)
library(rgeos)
library(plyr)
library(ggplot2)
library(reshape2)
library(maps)
library(mapproj)
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
################################################################################
# Read ESRI Shape File into SpatialPolygonsDataFrame
# rgdal readOGR() method gets projection information.
EcoRegion <- readOGR("/home/rorshach/Desktop/Ecomaps/Level3_state", "us_eco_l3_state_boundaries")
EcoRegion2<- spTransform(EcoRegion, CRS("+proj=longlat +datum=WGS84"))
################################################################################
EcoLegend=colsplit(EcoRegion2@data$L3_KEY," ",c("Code","name"))
EcoLegend2=ddply(EcoLegend, ~ name, summarise, Code= min(Code))
################################################################################
EcoRegion2@data$id = rownames(EcoRegion2@data)
EcoRegion2.points = fortify(EcoRegion2, region="id")
EcoRegion2.df = join(EcoRegion2.points, EcoRegion2@data, by="id")
###############################################################################
coords=as.data.frame(cbind(newDAPClist[,colnames(newDAPClist)=="Longitude"],newDAPClist[,colnames(newDAPClist)=="Latitude"]))
for (i in 1:nrow(coords)){
  if(is.na(coords[i,1])) newDAPClist$mapcode[i]="Unknown"
     else newDAPClist$mapcode[i]=as.character(over(SpatialPoints(coords[i,],proj4string=CRS(proj4string(EcoRegion2))), EcoRegion2)$L3_KEY)
    print(i)
    }
coords2=as.data.frame(cbind(newDAPClist[,colnames(newDAPClist)=="Longitude"],newDAPClist[,colnames(newDAPClist)=="Latitude"]))
for (i in 1:nrow(coords2)){
  if(is.na(coords2[i,1])) newDAPClist$mapcode[i]="Unknown"
  else newDAPClist$mapcode2[i]=as.character(over(SpatialPoints(coords2[i,],proj4string=CRS(proj4string(EcoRegion2))), EcoRegion2)$US_L3CODE)
  print(i)
}


temp=EcoLegend2[attributes,]
# identify some interesting attributes
attributes <- c("US_L3CODE","US_L3NAME","L3_KEY")
# subset the full dataset extracting only the desired attributes
SubsetEcoRegion <- EcoRegion2[,attributes]

# assign these attributes of interest to more descriptive names
#names(dataProjectedSubset) <- c("number", "area", "name")

# create a data.frame name (potentially different from layerName)
dataName <- "USECO"

# reproject the data onto a "longlat" projection and assign it to the new name
assign(dataName,spTransform(SubsetEcoRegion, CRS("+proj=longlat")))

# save the data
save(list=c(dataName),file=paste(localDir,"WAWRIAs.RData",sep="/"))

# inspect the watershed names
USECO$L3_KEY

PCGEcoNo=c(as.character(my3ddataset$mapcode))
PCGEcoSet=USECO[USECO$L3_KEY%in%PCGEcoNo,]
plot(PCGEcoSet) 



# create a "mask" identifying the biggest area
regionMask <- which(watershedSubsetDF$area == max(watershedSubsetDF$area))

#creat a subset
selectregion <- watershedSubsetDF[biggestAreaMask,]
################################################################################
EcoRegion2@data$subset =EcoRegion2@data$US_L3CODE
EcoRegion2[!EcoRegion2$L3_KEY%in%PCGEcoNo,colnames(EcoRegion2@data)=="subset"]=as.factor(c("1"))
EcoRegion2@data$subset

cluster = data.frame(id=EcoRegion2@data$id,
                     EcoRegion2@data$US_L3NAME,
                     cluster=EcoRegion2@data$subset)
EcoRegionplot=merge(EcoRegion2.df,cluster,by="id")
nlevels(EcoRegionplot$cluster)


ggplot(EcoRegionplot)+
    geom_polygon(aes(long,lat,group=group),fill = "white")+
  geom_path(aes(x=long,y=lat,group=group))+
  coord_equal()+
 theme_bw()+
geom_jitter(data=my3ddataset, aes(x=Longitude, y=Latitude,shape=clade, color=mapcode2,group=mapcode2))+
  scale_colour_manual(values = c("1"="white",
                             "27"="yellowgreen",
                             "28"="magenta4",
                             "40"="blue",
                             "43"="darkslategrey",
                             "46"="lightseagreen",
                             "47"= "darkgoldenrod3",
                             "48"="green3","cyan",
                             "51"="burlywood4",
                             "52"="darkred",
                             "53"="orangered",
                             "54"="purple",
                             "55"="ivory4",
                             "59"="plum4",
                             "72"="springgreen2",
                             "82"="hotpink4",
                             "84"="lightskyblue4"
                        )
)+
scale_shape_manual(values=c(0,1,2,15,16,17))


all_states <- map_data("state")
Qvalue=read.csv(file.choose())
mapdata=cbind(newDAPClist,Qvalue[,3:4])
mapdata[,8]=as.character(mapdata[,8])
mapdata[,9]=as.character(mapdata[,9])
mapdata[50,8]='KS'
mapdata[51,8]='KS'
mapdata[52,8]='KS'
mapdata[53,8]='KS'

mapdata[96,8]='Red River'
mapdata[96,9]='Red River'

gall <- subset(mapdata, State == "Red River")
gall_1 <- subset(mapdata, State2 == "NY2")
gall_2 <- subset(mapdata, State2 == "NY1")
gallnew=rbind(gall,gall_1,gall_2)
mapdata=mapdata[-c(1,2,96),]

png(file="ECOdata1.png",width=12.57703333,height=8,units="in",res=700)
p=ggplot(all_states)+geom_polygon(aes(long,lat,group=group),color="black", fill="white", size=0.25)+
  geom_path(aes(x=long,y=lat,group=group),color="black")+
  coord_equal()+
  theme_bw()+
  theme(axis.text= element_text(face = "bold", color = "black", size = 10))+
  theme(legend.justification=c(0,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_jitter(data=mapdata, size=6,aes(x=Longitude, y=Latitude,fill=deme,shape=as.factor(Ploidy)))+
  coord_map("albers", lat0=30, lat1=40)+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  #scale_fill_gradient(name = "%Membership",low="blue",high="red",limit=c(0,100))+
  scale_fill_manual(name = "Demes",values =c("blue","red"), labels=c("Deme1", "Deme2"))+
  guides(fill = guide_legend(override.aes = list(shape=18,color =c("blue","red"),size=8),order=1))+
  guides(shape = guide_legend(override.aes = list(size=5,stroke=1.5),order=2))+
  geom_point(data=gall,x=gall$Longitude,y=gall$Latitude, aes(colour="Red River"), shape=c(18), size=6,stroke = 2.3)+
  geom_point(data=gall_2[1,], x=gall_2$Longitude,y=gall_2$Latitude,aes(colour="KST"), shape=c(7), size=6,stroke = 1) +
  geom_point(data=gall_1, x=gall_1$Longitude,y=gall_1$Latitude,aes(colour="STP"), shape=c(9), size=6,stroke = 1)
#geom_label(data=gall,x=gall$Longitude,y=gall$Latitude,aes(label="Red River"),hjust=0,vjust=0,label.padding=unit(0.1, "lines"))+
#geom_label(data=gall_2[1,],x=gall_2$Longitude,y=gall_2$Latitude,aes(label="NY1"),hjust=0,vjust=0,label.padding=unit(0.1, "lines"))+
#geom_label(data=gall_1,x=gall_1$Longitude,y=gall_1$Latitude,aes(label="NY2"),hjust=0,vjust=0,label.padding=unit(0.1, "lines"))+
p=p+xlab("Longitude")+
  ylab("Latitude")+
  scale_colour_manual(name="Released varieties",values =c("black","black","orange"))+
  guides(colour = guide_legend(override.aes = list(shape=c(7,9,18),size=c(6,5,8),color =c("black","black","orange"))),order=3)
p
dev.off()


png(file="ECOdata2.png",width=12.57703333,height=8,units="in",res=700)
p=ggplot(all_states)+geom_polygon(aes(long,lat,group=group),color="black", fill="white", size=0.25)+
  geom_path(aes(x=long,y=lat,group=group),color="black")+
  coord_equal()+
  theme_bw()+
  theme(axis.text= element_text(face = "bold", color = "black", size = 10))+
  theme(legend.justification=c(0,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_jitter(data=mapdata, size=5,aes(x=Longitude, y=Latitude,fill=100*fast2,shape=as.factor(Ploidy)))+
  coord_map("albers", lat0=30, lat1=40)+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  scale_fill_gradient(name = "%Membership",low="blue",high="red",limit=c(0,100))+
  #scale_fill_manual(name = "Demes",values =c("Purple","orange"), labels=c("Deme1", "Deme2"))+
  #guides(fill = guide_legend(override.aes = list(shape=18,color =c("Purple","orange")),order=1))+
  guides(fill = guide_colorbar(label.theme = element_text(angle=0,size=12,
                                                          face="bold"),order=1))+
  guides(shape = guide_legend(override.aes = list(size=5,stroke=1.5),order=2))+
  geom_point(data=gallnew,aes(x=gallnew$Longitude,y=gallnew$Latitude,color=State2),fill=c("orange","black","black"),shape=c(7,9,18),size=c(4,5,6),stroke = 1)
  
  #geom_label(data=gall,x=gall$Longitude,y=gall$Latitude,aes(label="Red River"),hjust=0,vjust=0,label.padding=unit(0.1, "lines"))+
  #geom_label(data=gall_2[1,],x=gall_2$Longitude,y=gall_2$Latitude,aes(label="NY1"),hjust=0,vjust=0,label.padding=unit(0.1, "lines"))+
  #geom_label(data=gall_1,x=gall_1$Longitude,y=gall_1$Latitude,aes(label="NY2"),hjust=0,vjust=0,label.padding=unit(0.1, "lines"))+
p=p+xlab("Longitude")+
  ylab("Latitude")+
   scale_colour_manual(name="Released varieties",values =c("black","black","orange"))+
  guides(colour = guide_legend(override.aes = list(shape=c(7,9,18),size=c(6,5,8),color =c("black","black","orange"))),order=3)
p
 dev.off()

