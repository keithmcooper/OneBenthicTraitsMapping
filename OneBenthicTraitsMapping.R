################################################################################
####          CREATION OF RESPONSE/EFFECTS TRAIT RASTER LAYERS              ####
################################################################################

## This script relates to the work in Bolam, S.G., Cooper, K.M., Downie, A-L.,
# 2023. Mapping marine benthic biological traits to facilitate future
# sustainable development. Ecological Applications (in prep).

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/')

#_______________________________________________________________________________
## Load packages
library(ggplot2)
library(rgdal)
library(maptools)
library(plyr) 
library(pool)
library(DBI)
library (RPostgres)
library(dplyr)
library(lazyraster)
#_______________________________________________________________________________
#### CREATE A CONNECTION TO OneBenthic LIVE ####
Sys.setenv(R_CONFIG_ACTIVE = "one_benthic")

dw <- config::get()

pool <- dbPool(drv = dbDriver(dw$driver),
               dbname = dw$database,
               host = dw$server,
               port =  dw$port,
               user = dw$uid,
               password = dw$pwd)
#_______________________________________________________________________________
#### RETRIEVE DATA FROM ONEBENTHIC DB (CEFAS INTERNAL USE ONLY) #### 

## If you are not a member of Cefas staff, you should access the data from .csv (see below).

# CHOOSE EITHER:

## 1.SQL select query RESPONSE TRAITS
data = dbGetQuery(pool,
                  "SELECT
s.samplecode,
s.samplelat,
s.samplelong,
wtr.traits_modality,
SUM(wtr.fuzzyscore * sqrt(sqrt(ts.abund))) as traitscore,
SUM  (sqrt(sqrt(ts.abund))) as total_trans_abund

FROM samples.sample as s 
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 
INNER JOIN faunal_data.worrmstraits as wtr ON wtr.worrms_aphiaid = w.aphiaid

WHERE (wtr.traits_trait = 'morphology               ' OR wtr.traits_trait = 'eggdevelopment           ' OR wtr.traits_trait = 'livinghabit              ' OR wtr.traits_trait = 'sedimentposition         ' OR wtr.traits_trait = 'mobility                 ')
AND (s.gear_gearcode = 'MHN' OR s.gear_gearcode = 'DG' OR s.gear_gearcode = 'VV' OR s.gear_gearcode = 'SM' OR s.gear_gearcode = 'NIOZ' OR s.gear_gearcode = 'BC_0.1' OR s.gear_gearcode = 'C/VV'OR s.gear_gearcode = 'BC')
AND (s.treatment = 'R' or s.treatment IS NULL)
AND s.macrosieve= 1
AND w.include = TRUE
AND w.include_traits = TRUE
AND s.samplelat > 47.92938
AND s.id <= 44362
GROUP BY
s.samplecode,
wtr.traits_modality
ORDER by s.samplecode;")

# OR:

##2. SQL select query EFFECTS TRAITS
data = dbGetQuery(pool,
                  "SELECT
s.samplecode,
s.samplelat,
s.samplelong,
wtr.traits_modality,
SUM(wtr.fuzzyscore * sqrt(sqrt(ts.abund))) as traitscore,
SUM  (sqrt(sqrt(ts.abund))) as total_trans_abund

FROM samples.sample as s 
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 
INNER JOIN faunal_data.worrmstraits as wtr ON wtr.worrms_aphiaid = w.aphiaid

WHERE (wtr.traits_trait = 'maxsize                  ' OR wtr.traits_trait = 'longevity                ' OR wtr.traits_trait = 'larvaldevelopment        ' OR wtr.traits_trait = 'feedingmode              ' OR wtr.traits_trait = 'bioturbation             ')
AND (s.gear_gearcode = 'MHN' OR s.gear_gearcode = 'DG' OR s.gear_gearcode = 'VV' OR s.gear_gearcode = 'SM' OR s.gear_gearcode = 'NIOZ' OR s.gear_gearcode = 'BC_0.1' OR s.gear_gearcode = 'C/VV'OR s.gear_gearcode = 'BC')
AND (s.treatment = 'R' or s.treatment IS NULL)
AND s.macrosieve= 1
AND w.include = TRUE
AND w.include_traits = TRUE
AND s.samplelat > 47.92938
AND s.id <= 44362
GROUP BY
s.samplecode,
wtr.traits_modality
ORDER by s.samplecode;")

#_______________________________________________________________________________
#### RETRIEVE DATA FROM .csv files (see https://doi.org/10.14466/CefasDataHub.137) #### 

## Set working directory
setwd("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTraitsMapping")

# CHOOSE EITHER:

## 1.Load RESPONSE TRAITS dataset from .csv file
data=read.csv("DATA/response_public.csv", header=T,na.strings=c("NA", "-","?","<null>"),
                  stringsAsFactors=F,check.names=FALSE)

# OR:

## 2.Load EFFECTS TRAITS dataset from .csv file
data=read.csv("DATA/effects_public.csv", header=T,na.strings=c("NA", "-","?","<null>"),
              stringsAsFactors=F,check.names=FALSE)

## Drop 1st column
data1 <- data1[,2:7]
#_______________________________________________________________________________
## Check data
head(data)

## Check the number of samples
df_uniq <- unique(data$samplecode)
length(df_uniq) #31838

## Calculate traitscore as a percentage of total transformed abundance
data$tperc <- (data$traitscore/data$total_trans_abund)*100

# Drop unwanted columns
data <- data[,c(1,2,3,4,7)]
colnames(data)[5] <- "traitscore"

#_______________________________________________________________________________
#### PREPARE TRAITS DATA FOR CLUSTERING ####

## Examine col names of retreived data
head(data)# 1.samplecode, 2. samplelat, 3. samplelong, 4. traits_modality, 5. traitscore

## Load library
library(tidyr)

## Change data from long to wide format: data, key (headers), values
data2 <-spread(data,traits_modality, traitscore)
dim(data2) #31838    27
head(data2)
#_______________________________________________________________________________
#### 10. REMOVE 'REPLICATES' TO ADDRESS SPATIAL AUTOCORRELATION ####

## Load library
library(sp)

## Set coordinates
coordinates(data2) <- c("samplelong", "samplelat")

## Work out 50m distance in decimal degrees. 1 degree of latitude =111,000m.
50/111000# degrees for 50m #0.0004504505

## Set distance within which to remove replicates (zero)
zd <- zerodist(data2,zero = 0.0004504505)

## Drop replicates
data3 <- data2[-zd[,2], ]
dim(data3)# 18336    25

## Change class to df
data3 <- data.frame(data3)
class(data3)
names(data3)
head(data3)

## Drop col 'optional'
data4 <- data3[,1:27]
head(data4)

## Plot sample positions
plot(data3$samplelong,data3$samplelat)

#_______________________________________________________________________________
#### PREPARE DATA FOR CLUSTERING ####

## Select only the trait modality columns
data5=data4[,4:ncol(data4)]

## Check dimensions of df 'data5'
dim(data5) #18336    24

## Check df 'data5' is just the traits data 
names(data5)# it is 

## Change class of df data5 to a matrix 
data6=data.matrix(data5) 

## Create a df 'pos' for samplecode, samplelat, samplelong
pos=data4[,1:3] 

## Check names
names(pos)

## Make df
pos <- as.data.frame(pos)
class(pos)
#_______________________________________________________________________________
#### TRANSFORM DATA ####

## 4th root transformation applied in initial sql query therefore not required
#datat=data6^(0.25)
datat=data6
#_______________________________________________________________________________
#### ELBOW PLOT TO IDENTIFY APPROPRIATE NUMBER OF CLUSTER GROUPS ####
#see # http://www.r-statistics.com/2013/08/k-means-clustering-from-r-in-action/) 

wssplot <- function(datat, nc=40, seed=1234){ 
  wss <- (nrow(datat)-1)*sum(apply(datat,2,var)) 
  for (i in 2:nc){ 
    set.seed(seed) 
    wss[i] <- sum(kmeans(datat, algorithm="MacQueen",iter.max=100,nstart=25,centers=i)$withinss)} 
  plot(1:nc, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")} 

## Save plot to an image file (png of tiff) 
# Tosave as tiff use '.tiff' instead of '.png'

# CHOOSE EITHER:
png('OUTPUTS/RESPONSE_TRAITS/rt_elbow.png') #width = 15,height = 13,units = "cm", res = 800, pointsize = 12

#OR: 
png('OUTPUTS/EFFECTS_TRAITS/et_elbow.png') #width = 15,height = 13,units = "cm", res = 800, pointsize = 12

## Set text size 
par(cex=1) #0.8

## Generate plot 
wssplot(datat)

## Add line to show position of elbow 
abline(v=6,col="red")

dev.off()
#_______________________________________________________________________________
#### CLUSTERING OF TRAITS DATA ####

set.seed(1234)

## Perform clustering
results=kmeans(datat,6,algorithm="MacQueen",iter.max=100,nstart=25)
results

#_______________________________________________________________________________
#### DENDROGRAM ####

## Calculate the absolute differences between cluster centres over all variables.
nclusters = 6
absdiff = matrix(0, nrow=nclusters, ncol=nclusters)
centers = results$centers
for (j in 1:nclusters) {
  for (k in 1:nclusters) {
    absdiff[j,k] = sum(abs(centers[j,] - centers[k,]))
  }
}
d=round(absdiff, 1)

## Find distance matrix
d1 <- dist(as.matrix(d))

## Apply hirarchical clustering
hc <- hclust(d1)

## Save plot of dendrogram as an image file(png or tiff)
# CHOOSE
#png("OUTPUTS/RESPONSE_TRAITS/rt_dendro.png", width = 15, height = 13, units = "cm", res = 800,pointsize = 12) # RESPONSE TRAITS OR
png("OUTPUTS/EFFECTS_TRAITS/et_dendro.png", width = 15, height = 13, units = "cm", res = 800,pointsize = 12) # EFFECTS TRAITS

## Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(hc)

## Reduce label size of cluster labels 1-5
#par(cex=0.8)
par(cex=0.8)

## Plot dendrogram
plot(hcd, ylab = "Height",leaflab = "none")
#plot(hcd, ylab = "Height")# with labels

## Add circle symbols below each leaf of the dendrogram. Define colours.Useful sites for colours: https://colorbrewer2.org/#type=diverging&scheme=PRGn&n=6 or 
#col.circle=c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59') # RESPONE TRAITS https://mycolor.space/?hex=%23BAC25E&sub=1 ##F6906C
col.circle=c('#FF0000', '#c7eae5','#5ab4ac', '#d8b365', '#8c510a', '#f6e8c3') # EFFECTS TRAITS https://mycolor.space/?hex=%23BAC25E&sub=1

## Add symbols
symbols(1:6, rep(0, 6), circles=rep(1, 6), add=TRUE, inches=.08,fg=col.circle,
        bg=col.circle, xpd=TRUE)#-10

## Add cluster group labels
#axis(1,at=seq(1,6,by=1),labels=c("2","1","3","4","5","6"),pos=-5,cex.axis=0.8,lty = 0,cex.axis=1)# RESPONSE
axis(1,at=seq(1,6,by=1),labels=c("6","3","1","2","4","5"),pos=-5,cex.axis=0.8,lty = 0,cex.axis=1)# EFFECTS

dev.off()
#_______________________________________________________________________________
#### PREPARE CLUSTER RESULTS FOR PLOTTING ####

library(ggplot2)
library(mapproj)

## Join clustering results to sample positions
samclus=cbind(pos,results$cluster)
head(samclus)

##Rename cluster results column
colnames(samclus)[4] <- 'cluster'
str(samclus)
samclus$cluster <- as.factor(samclus$cluster)

## Plot
ggplot()+ 
  geom_point(data=samclus,aes(samplelong,samplelat,col=cluster), size=0.6,show.legend = TRUE)+
  # CHOOSE EITHER
  scale_colour_manual(values = c('#4575b4', '#00E600', '#91bfdb', '#e0f3f8','#fee090','#fc8d59'),name="Cluster")+# RESPONSE TRAITS OR
  #scale_colour_manual(values = c('#5ab4ac', '#d8b365', '#c7eae5', '#8c510a','#f6e8c3', '#FF0000'),name="Cluster")+# EFFECTS TRAITS
  #geom_polygon(data=euDF2, aes(x=long, y=lat, group=group),fill="white",colour="black", size=0.05)+ 
  coord_map(xlim = c(-10, 10),ylim = c(48, 61))+
  guides(colour = guide_legend(override.aes = list(size=3))) # Change size of legend dots
#_______________________________________________________________________________
#### PLOT TRAIT CLUSTERS USING LEAFLET ####

library(leaflet)
library(leafem)

# Define cluster colours
pal <- colorFactor(
  palette = c('#4575b4', '#00E600', '#91bfdb', '#e0f3f8','#fee090','#fc8d59'),domain = samclus$cluster)# RESPONSE TRAITS
#palette = c('#5ab4ac', '#d8b365', '#c7eae5', '#8c510a','#f6e8c3', '#FF0000'),domain = samclus$cluster)# EFFECTS TRAITS

## Map
leaflet() %>%
  addProviderTiles(providers$Esri.OceanBasemap,options = providerTileOptions(noWrap = TRUE))%>%
  #addPolygons(data=licence, weight = 1,fillColor="white",fillOpacity = 0)%>%
  #addPolygons(data=mcz, weight = 1,color="orange",fillColor="orange",fillOpacity = 0)%>%
  addCircleMarkers(data=samclus,~as.numeric(samplelong), ~as.numeric(samplelat), popup = ~as.character(cluster),radius = 3,stroke = F, color = "black",weight = 1,fill = TRUE, fillColor =~pal(cluster),fillOpacity = 1)%>%
  setView(-3,54.6,zoom=5)%>%
  addMouseCoordinates()
#_______________________________________________________________________________
#### CLUSTER PLOTS ####
##https://www.guru99.com/r-k-means-clustering.html

## Get cluster centres
center <-results$centers
head(center)
names(as.data.frame(center))

library(tidyr)

## Add cluster group
cluster <- c(1:6)

## Create dataset with the cluster number
center_df <- data.frame(cluster, center)

# Reshape the data
# Choose either
center_reshape <- gather(center_df, features, values, edAsexual.Budding: spSurface)# RESPONSE
#center_reshape <- gather(center_df, features, values, bDiffusive.mixing: sr21.100)# EFFECTS

head(center_reshape)

library(RColorBrewer)

## Create the palette
hm.palette <-colorRampPalette(rev(brewer.pal(10, 'RdYlGn')),space='Lab')

## Change order of traits to match gt table
# RESPONSE TRAITS
#center_reshape$features <- factor(center_reshape$features, levels = rev(c("mSoft","mTunic","mExoskeleton..chitin.calcium.carbonate.","mCrustose","mCushion","mStalked","edAsexual.Budding","edSexual.shed.eggs..pelagic","edSexual.shed.eggs..benthic","edSexual.brood.eggs","lhTube.dwelling","lhBurrow.dwelling","lhFree.living","lhCrevice.hole.under.stone","lhEpi.endo.biotic","lhAttached.to.substratum","spSurface","spInfauna..0.5cm","spInfauna..6.10cm","spInfauna...10cm","mobSessile","mobCrawl.creep.climb","mobBurrower","mobSwim")))
# EFFECTS TRAITS
center_reshape$features <- factor(center_reshape$features, levels = rev(c("sr.10","sr11.20","sr21.100","sr101.200","sr201.500","sr.500","l.1","l1.2..should.be..1...3.","l3.10","l.10","ldPlanktotrophic","ldLecithotrophic","ldDirect","fSuspension","fSurface.Deposit","fSubsurface.deposit","fScavenger.Opportunist","fPredator","fParasite","bDiffusive.mixing","bSurface.deposition","bUpward.Conveyor","bDownwards.conveyer","bNone")))


## Plot the heat map
p <- ggplot(data = center_reshape, aes(x = features, y = cluster, fill = values)) +
  #scale_y_continuous(breaks = seq(1, 6, by = 1)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()#+
#theme(axis.text.x = element_text(angle = 90))

#ggsave(p2, file="OUTPUTS/RESPONSE_TRAITS/rt_cluster_plot2.png" , width=10, height=8)# RESPONSE TRAITS
#p2 <- p+coord_flip()
#png('OUTPUTS/RESPONSE_TRAITS/rt_cluster_plot.png')
png('OUTPUTS/EFFECTS_TRAITS/et_cluster_plot.png')
p+coord_flip()
dev.off()
#_______________________________________________________________________________
#### PREPARE DATA FOR ADDING TO GT TABLE

## Require df traits, traits_modality, 1, 2, 3, 4, 5
head(center_reshape)

## Long to wide format
library(tidyr)
centers4gt <- center_reshape %>% spread(cluster, values, fill = NA, convert = FALSE)
class(centers4gt)

## Add new col for traits
centers4gt$traits <- centers4gt$features
head(centers4gt)

## Reorder columns
centers4gt2 <- centers4gt[,c(8,1:7)]

## Update column names
colnames(centers4gt2) <- c("traits","traits_modality","1","2","3","4","5","6")
head(centers4gt2)
unique(centers4gt2$traits)
unique(centers4gt2$traits_modality)
str(centers4gt2)

## Change from factor to chr (otherwise you can change names in next step)
centers4gt2$traits <- as.character(centers4gt2$traits)
centers4gt2$traits_modality <- as.character(centers4gt2$traits_modality)

## Update names of traits_modality modalities
# Morphology
centers4gt2$traits_modality[centers4gt2$traits_modality == "mCrustose"] <- "Crustose"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mCushion"] <- "Cushion"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mExoskeleton..chitin.calcium.carbonate."] <- "Exoskeleton"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mSoft" ] <- "Soft"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mStalked"] <- "Stalked"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mTunic"] <- "Tunic"

# Egg development
centers4gt2$traits_modality[centers4gt2$traits_modality == "edAsexual.Budding"] <- "Asexual"
centers4gt2$traits_modality[centers4gt2$traits_modality == "edSexual.brood.eggs"] <- "Brood"
centers4gt2$traits_modality[centers4gt2$traits_modality == "edSexual.shed.eggs..benthic"] <- "Benthic"
centers4gt2$traits_modality[centers4gt2$traits_modality == "edSexual.shed.eggs..pelagic"] <- "Pelagic"

# Living habit
centers4gt2$traits_modality[centers4gt2$traits_modality == "lhAttached.to.substratum"] <- "Attached"
centers4gt2$traits_modality[centers4gt2$traits_modality == "lhBurrow.dwelling"] <- "Burrow"
centers4gt2$traits_modality[centers4gt2$traits_modality == "lhCrevice.hole.under.stone"] <- "Crevice"
centers4gt2$traits_modality[centers4gt2$traits_modality == "lhEpi.endo.biotic"] <- "Epi/endo/zoic/phytic"
centers4gt2$traits_modality[centers4gt2$traits_modality == "lhFree.living"] <- "Free-living"
centers4gt2$traits_modality[centers4gt2$traits_modality == "lhTube.dwelling"] <- "Tube"

# Sediment position
centers4gt2$traits_modality[centers4gt2$traits_modality == "spInfauna...10cm"] <- ">10cm"
centers4gt2$traits_modality[centers4gt2$traits_modality == "spInfauna..0.5cm"] <- "0-5cm"
centers4gt2$traits_modality[centers4gt2$traits_modality == "spInfauna..6.10cm"] <- "6-10cm"
centers4gt2$traits_modality[centers4gt2$traits_modality == "spSurface"] <- "Surface"

# Mobility
centers4gt2$traits_modality[centers4gt2$traits_modality == "mobBurrower"] <- "Burrower"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mobCrawl.creep.climb"] <- "Crawl/creep/climb"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mobSessile"] <- "Sessile"
centers4gt2$traits_modality[centers4gt2$traits_modality == "mobSwim"] <- "Swim"

# Body size
centers4gt2$traits_modality[centers4gt2$traits_modality == "sr.10"] <- "<10"
centers4gt2$traits_modality[centers4gt2$traits_modality == "sr.500"] <- "500"
centers4gt2$traits_modality[centers4gt2$traits_modality == "sr101.200"] <- "101-200"
centers4gt2$traits_modality[centers4gt2$traits_modality == "sr11.20"] <- "11-20"
centers4gt2$traits_modality[centers4gt2$traits_modality == "sr201.500"] <- "201-500"
centers4gt2$traits_modality[centers4gt2$traits_modality == "sr21.100"] <- "21-100"

# Longevity
centers4gt2$traits_modality[centers4gt2$traits_modality == "l.1"] <- "<1"
centers4gt2$traits_modality[centers4gt2$traits_modality == "l.10"] <- ">10"
centers4gt2$traits_modality[centers4gt2$traits_modality == "l1.2..should.be..1...3."] <- "1-2"
centers4gt2$traits_modality[centers4gt2$traits_modality == "l3.10"] <- "3-10"

# Larval development
centers4gt2$traits_modality[centers4gt2$traits_modality == "ldDirect"] <- "Benthic (Direct)"
centers4gt2$traits_modality[centers4gt2$traits_modality == "ldLecithotrophic"] <- "Lecithotrophic"
centers4gt2$traits_modality[centers4gt2$traits_modality == "ldPlanktotrophic"] <- "Planktotrophic"

# Feeding mode
centers4gt2$traits_modality[centers4gt2$traits_modality == "fParasite"] <- "Parasite"
centers4gt2$traits_modality[centers4gt2$traits_modality == "fPredator"] <- "Predator"
centers4gt2$traits_modality[centers4gt2$traits_modality == "fScavenger.Opportunist"] <- "Scavenger/Opportunist"
centers4gt2$traits_modality[centers4gt2$traits_modality == "fSubsurface.deposit"] <- "Sub-surface Deposit"
centers4gt2$traits_modality[centers4gt2$traits_modality == "fSurface.Deposit"] <- "Surface Deposit"
centers4gt2$traits_modality[centers4gt2$traits_modality == "fSuspension"] <- "Suspension"

# Bioturbation
centers4gt2$traits_modality[centers4gt2$traits_modality == "bDiffusive.mixing"] <- "Diffusive Mixing"
centers4gt2$traits_modality[centers4gt2$traits_modality == "bDownwards.conveyer"] <- "Downward Conveyer"
centers4gt2$traits_modality[centers4gt2$traits_modality == "bNone"] <- "None"
centers4gt2$traits_modality[centers4gt2$traits_modality == "bSurface.deposition"] <- "Surface Deposition"
centers4gt2$traits_modality[centers4gt2$traits_modality == "bUpward.Conveyor"] <- "Upward Conveyor"

## update traits names
# Morphology
centers4gt2$traits[centers4gt2$traits == "mCrustose"] <- "Morphology"
centers4gt2$traits[centers4gt2$traits == "mCushion"] <- "Morphology"
centers4gt2$traits[centers4gt2$traits == "mExoskeleton..chitin.calcium.carbonate."] <- "Morphology"
centers4gt2$traits[centers4gt2$traits == "mSoft"] <- "Morphology"
centers4gt2$traits[centers4gt2$traits == "mStalked"] <- "Morphology"
centers4gt2$traits[centers4gt2$traits == "mTunic"] <- "Morphology"

# Egg development
centers4gt2$traits[centers4gt2$traits == "edAsexual.Budding"] <- "Egg Development"
centers4gt2$traits[centers4gt2$traits == "edSexual.brood.eggs"] <- "Egg Development"
centers4gt2$traits[centers4gt2$traits == "edSexual.shed.eggs..benthic"] <- "Egg Development"
centers4gt2$traits[centers4gt2$traits == "edSexual.shed.eggs..pelagic"] <- "Egg Development"

# Living habit
centers4gt2$traits[centers4gt2$traits == "lhAttached.to.substratum"] <- "Living Habit"
centers4gt2$traits[centers4gt2$traits == "lhBurrow.dwelling"] <- "Living Habit"
centers4gt2$traits[centers4gt2$traits == "lhCrevice.hole.under.stone"] <- "Living Habit"
centers4gt2$traits[centers4gt2$traits == "lhEpi.endo.biotic"] <- "Living Habit"
centers4gt2$traits[centers4gt2$traits == "lhFree.living"] <- "Living Habit"
centers4gt2$traits[centers4gt2$traits == "lhTube.dwelling"] <- "Living Habit"

# Sediment position
centers4gt2$traits[centers4gt2$traits == "spInfauna...10cm"] <- "Sediment Position"
centers4gt2$traits[centers4gt2$traits == "spInfauna..0.5cm"] <- "Sediment Position"
centers4gt2$traits[centers4gt2$traits == "spInfauna..6.10cm"] <- "Sediment Position"
centers4gt2$traits[centers4gt2$traits == "spSurface"] <- "Sediment Position"

# Mobility
centers4gt2$traits[centers4gt2$traits == "mobBurrower"] <- "Mobility"
centers4gt2$traits[centers4gt2$traits == "mobCrawl.creep.climb"] <- "Mobility"
centers4gt2$traits[centers4gt2$traits == "mobSessile"] <- "Mobility"
centers4gt2$traits[centers4gt2$traits == "mobSwim"] <- "Mobility"

# Body size
centers4gt2$traits[centers4gt2$traits == "sr.10"] <- "Maxsize"
centers4gt2$traits[centers4gt2$traits == "sr.500"] <- "Maxsize"
centers4gt2$traits[centers4gt2$traits == "sr101.200"] <- "Maxsize"
centers4gt2$traits[centers4gt2$traits == "sr11.20"] <- "Maxsize"
centers4gt2$traits[centers4gt2$traits == "sr201.500"] <- "Maxsize"
centers4gt2$traits[centers4gt2$traits == "sr21.100"] <- "Maxsize"

# Longevity
centers4gt2$traits[centers4gt2$traits == "l.1"] <- "Longevity"
centers4gt2$traits[centers4gt2$traits == "l.10"] <- "Longevity"
centers4gt2$traits[centers4gt2$traits == "l1.2..should.be..1...3."] <- "Longevity"
centers4gt2$traits[centers4gt2$traits == "l3.10"] <- "Longevity"

# Larval development
centers4gt2$traits[centers4gt2$traits == "ldDirect"] <- "Larval Development"
centers4gt2$traits[centers4gt2$traits == "ldLecithotrophic"] <- "Larval Development"
centers4gt2$traits[centers4gt2$traits == "ldPlanktotrophic"] <- "Larval Development"

# Feeding mode
centers4gt2$traits[centers4gt2$traits == "fParasite"] <- "Feeding Mode"
centers4gt2$traits[centers4gt2$traits == "fPredator"] <- "Feeding Mode"
centers4gt2$traits[centers4gt2$traits == "fScavenger.Opportunist"] <- "Feeding Mode"
centers4gt2$traits[centers4gt2$traits == "fSubsurface.deposit"] <- "Feeding Mode"
centers4gt2$traits[centers4gt2$traits == "fSurface.Deposit"] <- "Feeding Mode"
centers4gt2$traits[centers4gt2$traits == "fSuspension"] <- "Feeding Mode"

# Bioturbation
centers4gt2$traits[centers4gt2$traits == "bDiffusive.mixing"] <- "Bioturbation"
centers4gt2$traits[centers4gt2$traits == "bDownwards.conveyer"] <- "Bioturbation"
centers4gt2$traits[centers4gt2$traits == "bNone"] <- "Bioturbation"
centers4gt2$traits[centers4gt2$traits == "bSurface.deposition"] <- "Bioturbation"
centers4gt2$traits[centers4gt2$traits == "bUpward.Conveyor"] <- "Bioturbation"

## Check names changes successfully
head(centers4gt2)
unique(centers4gt2$traits)
unique(centers4gt2$traits_modality)
#_______________________________________________________________________________
#### CHANGE ORDER OF MODALITIES ####

## Define order or rows
desiredOrder <- c(
  "Soft",
  "Tunic",
  "Exoskeleton",
  "Crustose",
  "Cushion",
  "Stalked",
  "Asexual",
  "Pelagic",
  "Benthic",
  "Brood",
  "Tube",
  "Burrow",
  "Free-living",
  "Crevice",
  "Epi/endo/zoic/phytic",
  "Attached",
  "Surface",
  "0-5cm",
  "6-10cm",
  ">10cm",
  "Sessile",
  "Crawl/creep/climb",
  "Burrower",
  "Swim", 
  "<10",
  "11-20",
  "21-100",
  "101-200",              
  "201-500",
  "500",
  "<1",
  "1-2",
  "3-10",
  ">10",
  "Planktotrophic",
  "Lecithotrophic",
  "Benthic (Direct)",
  "Suspension",
  "Surface Deposit",
  "Sub-surface Deposit",
  "Scavenger/Opportunist",
  "Predator",
  "Parasite",
  "Diffusive Mixing",
  "Surface Deposition",
  "Upward Conveyor",
  "Downward Conveyer",
  "None"
)

## Update order based on above df
traits_mean_perc5 <- centers4gt2[match(desiredOrder, centers4gt2$traits_modality), ] 
str(traits_mean_perc5)
head(traits_mean_perc5)

## Remove NA rows
traits_mean_perc5 <-traits_mean_perc5[complete.cases(traits_mean_perc5), ]

#_______________________________________________________________________________
## Create table RESPONSE TRAITS
library(gt)
chara_traits <- traits_mean_perc5%>%gt()%>%
  fmt_number(columns = vars('1','2','3','4','5','6'), decimals = 1)%>%
  tab_spanner(
    label = md("**Cluster means**"),
    columns = vars('2','1','3','4','5','6')
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('1'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('2'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('3'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('4'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  ## Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('5'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  ## Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('6'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    ))%>%
  tab_style(
    style = list(
      cell_fill(color = "#fc8d59")
    ),
    location = list(
      cells_column_labels(columns = vars('6'))
    ))%>%
  tab_style(
    style = list(
      cell_fill(color = "#91bfdb")
    ),
    location = list(
      cells_column_labels(columns = vars('3'))
    ))%>% 
  tab_style(
    style = list(
      cell_fill(color = "#e0f3f8")
    ),
    location = list(
      cells_column_labels(columns = vars('4'))
    ))%>%
  # colour col header
  tab_style(
    style = list(
      cell_fill(color = "#fee090")#,
      #cell_text(color = "white")
    ),
    location = list(
      cells_column_labels(columns = vars('5'))
    ))%>%
  tab_style(
    style = list(
      cell_fill(color = "#4575b4")
    ),
    location = list(
      cells_column_labels(columns = vars('1'))
    ))%>%
  # colour col header Clus2
  tab_style(
    style = list(
      cell_fill(color = "#00E600")
    ),
    location = list(
      cells_column_labels(columns = vars('2'))
    ))%>% 
  # colour col header
  cols_label(
    traits_modality = md("**Modality**"),
    traits = md("**Trait**"),
    '1' = md("**1**"),
    '2' = md("**2**"),
    '3' = md("**3**"),
    '4' = md("**4**"),
    '5' = md("**5**"),
    '6' = md("**6**"))

## View table
chara_traits

## Save Response traits table
chara_traits %>%
  gtsave(
    "OUTPUTS/RESPONSE_TRAITS/rt_characterising_traits_tablev1.png")

#_______________________________________________________________________________
## Create table EFFECTS TRAITS
library(gt)
chara_traits <- traits_mean_perc5%>%gt()%>%
  fmt_number(columns = vars('1','2','3','4','5','6'), decimals = 1)%>%
  tab_spanner(
    label = md("**Cluster means**"),
    columns = vars('6','3','1','2','4','5')
  )%>%
  ## Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('6'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    ))%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('4'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('1'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('2'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  # Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('3'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  
  ## Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
  data_color( # Update cell colors...
    columns = vars('5'), # ...for dose column 
    colors = scales::col_numeric( # <- bc it's numeric
      palette = c(
        "white","#737373"), # A color scheme (gradient)
      domain = c(0,77) # Column scale endpoints
    )
  )%>%
  
  tab_style(
    style = list(
      cell_fill(color = '#5ab4ac')
    ),
    location = list(
      cells_column_labels(columns = vars('1'))
    ))%>%
  # colour col header Clus2
  tab_style(
    style = list(
      cell_fill(color ='#d8b365' )
    ),
    location = list(
      cells_column_labels(columns = vars('2'))
    ))%>% 
  
  tab_style(
    style = list(
      cell_fill(color = '#c7eae5')
    ),
    location = list(
      cells_column_labels(columns = vars('3'))
    ))%>% 
  
  tab_style(
    style = list(
      cell_fill(color = '#8c510a')
    ),
    location = list(
      cells_column_labels(columns = vars('4'))
    ))%>%
  # colour col header
  tab_style(
    style = list(
      cell_fill(color = '#f6e8c3')#,
      #cell_text(color = "white")
    ),
    location = list(
      cells_column_labels(columns = vars('5'))
    ))%>%
  tab_style(
    style = list(
      cell_fill(color = '#FF0000')
    ),
    location = list(
      cells_column_labels(columns = vars('6'))
    ))%>%
  # colour col header
  cols_label(
    traits_modality = md("**Modality**"),
    traits = md("**Trait**"),
    '6' = md("**6**"),
    '4' = md("**4**"),
    '1' = md("**1**"),
    '2' = md("**2**"),
    '3' = md("**3**"),
    '5' = md("**5**"))

## View table
chara_traits

## Save Effects traits table
chara_traits %>%
  gtsave(
    "OUTPUTS/EFFECTS_TRAITS/et_characterising_traits_table.png")
#_______________________________________________________________________________
#### PCA ####
# see https://cran.r-project.org/web/packages/factoextra/readme/README.html

library(FactoMineR)
library(factoextra)

## Start with the data used for clustering in df datat
head(datat)
class(datat)# matrix array

## Create a copy of datat
datat2 <- as.data.frame(datat)

## Update name of columns in so they are suitable for use in biplot
names(datat2)
# RESPONSE TRAITS
#colnames(datat2) <- c("edAsexual","edBrood","edBenthic","edPelagic","lhAttached","lhBurrow","lhCrevice","lhEpi","lhFree","lhTube","mCrustose","mCushion","mExoskeleton","mobBurrower",
#                      "mobCrawl","mobSessile","mobSwim","mSoft","mStalked","mTunic","sp>10cm","sp0-5cm","sp6-10cm","spSurface"  )
colnames(datat2) <- c("Asexual","Brood","Benthic","Pelagic","Attached","Burrow","Crevice","Epi/endo/zoic/phytic","Free-living","Tube","Crustose","Cushion","Exoskeleton","Burrower","Crawl/creep/climb",
                      "Sessile","Swim","Soft","Stalked","Tunic",">10cm","0-5cm","6-10cm","Surface")# RESPONSE TRAITS
colnames(datat2) <- c("Diffusive Mixing","Downwards Conveyer","None","Surface Deposit","Upward Conveyor","Parasite","Predator",
                      "Scavenger/Opportunist","Subsurface Deposition","Surface Deposition","Suspension","<1",">10","1-2","3-10",
                      "Benthic (Direct)","Lecithotrophic","Planktotrophic","<10","500","101-200","11-20","201-500","21-100")# EFFECTS TRAITS
## Turn object datat2 back to a matrix
datat2 <- as.matrix(datat2)
class(datat2)
View(datat2)

## PCA
res.pca <- PCA(datat2,  graph = FALSE)

# Extract eigenvalues/variances
get_eig(res.pca)

## Scree plot
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

## Extract and visualize results for variables
# Extract the results for variables
var <- get_pca_var(res.pca)
var

# Coordinates of variables
head(var$coord)

# Contribution of variables
head(var$contrib)

# Graph of variables: default plot
fviz_pca_var(res.pca, col.var = "black")

#png(file="OUTPUTS/RESPONSE_TRAITS/rt_pca_no_samples.png", width=500, height=500)#RESPONSE TRAITS
png(file="OUTPUTS/EFFECTS_TRAITS/et_pca_no_samples.png", width=500, height=500)# EFFECTS TRAITS
fviz_pca_var(res.pca, col.var = "black")
dev.off()

## Biplot (RESPONSE TRAITS)
#devtools::install_github("vqv/ggbiplot")
library(ggbiplot)


myplot <- res.pca %>% 
  ggbiplot::ggbiplot(
    scale = 1, 
    alpha=0.3,
    varname.abbrev=FALSE,
    var.axes = F,
    varname.size = 3.7,
    varname.adjust = 1.5,
    groups = factor(results$cluster),
    ellipse = F)+
  scale_color_manual(values=c('#4575b4', '#00E600', '#91bfdb', '#e0f3f8','#fee090','#fc8d59'),name="Cluster")+
  geom_vline(xintercept = 0, linetype = 3,size=1)+
  geom_hline(yintercept = 0, linetype = 3,size=1)+
  ggpubr::theme_pubr(border = TRUE,
                     margin = TRUE,) + 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),legend.text=element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(size=5,alpha = 1)))+
  labs(x="PC1 (33.7%)",y="PC2 (12.9%)")

myplot

## Biplot (EFFECTS TRAITS)
#devtools::install_github("vqv/ggbiplot")
library(ggbiplot)

myplot <- res.pca %>% 
  ggbiplot::ggbiplot(
    scale = 1, 
    alpha=0.3,
    varname.abbrev=FALSE,
    var.axes = F,
    varname.size = 3.7,
    varname.adjust = 1.5,
    groups = factor(results$cluster),
    ellipse = F)+
  scale_color_manual(values=c('#5ab4ac', '#d8b365', '#c7eae5', '#8c510a','#f6e8c3', '#FF0000'),name="Cluster")+
  geom_vline(xintercept = 0, linetype = 3,size=1)+
  geom_hline(yintercept = 0, linetype = 3,size=1)+
  ggpubr::theme_pubr(border = TRUE,
                     margin = TRUE,) + 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),legend.text=element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(size=5,alpha = 1)))+
  labs(x="PC1 (18.8%)",y="PC2 (14.5%)")

## Bring arrows to front 
myplot$layers <- c(myplot$layers, myplot$layers[[1]])
myplot

## Save plot
ggsave(myplot, file="OUTPUTS/RESPONSE_TRAITS/ResponseTraitsPCA_no_labs.png" , width=10, height=8)# RESPONSE TRAITS
ggsave(myplot, file="OUTPUTS/EFFECTS_TRAITS/EffectsTraitsPCA_no_labs.png" , width=10, height=8)# EFFECTS TRAITS
#_______________________________________________________________________________
#### PCA CONTRIBUTIONS TABLE ####

## Start with PCA coordinates
head(var$coord)
class(var$coord)

## Change from matrix to df
coord <- data.frame(Modality = row.names(var$coord), var$coord)

## Drop rownames
rownames(coord) <- NULL
coord

## Add col for trait
coord$Trait <- coord$Modality

## Now replace values in col trait with correct string
coord$Trait [coord$Trait  == "Soft"] <- "Morphology"
coord$Trait [coord$Trait  == "Tunic"] <- "Morphology"
coord$Trait [coord$Trait  == "Exoskeleton"] <- "Morphology"
coord$Trait [coord$Trait  == "Crustose"] <- "Morphology"
coord$Trait [coord$Trait  == "Cushion"] <- "Morphology"
coord$Trait [coord$Trait  == "Stalked"] <- "Morphology"

coord$Trait [coord$Trait  == "Asexual"] <- "Egg Development"
coord$Trait [coord$Trait  == "Pelagic"] <- "Egg Development"
coord$Trait [coord$Trait  == "Benthic"] <- "Egg Development"
coord$Trait [coord$Trait  == "Brood"] <- "Egg Development"

coord$Trait [coord$Trait  == "Tube"] <- "Living Habit"
coord$Trait [coord$Trait  == "Burrow"] <- "Living Habit"
coord$Trait [coord$Trait  == "Free-living"] <- "Living Habit"
coord$Trait [coord$Trait  == "Crevice"] <- "Living Habit"
coord$Trait [coord$Trait  == "Epi/endo/zoic/phytic"] <- "Living Habit"
coord$Trait [coord$Trait  == "Attached"] <- "Living Habit"

coord$Trait [coord$Trait  == "Surface"] <- "Sediment Position"
coord$Trait [coord$Trait  == "0-5cm"] <- "Sediment Position"
coord$Trait [coord$Trait  == "6-10cm"] <- "Sediment Position"
coord$Trait [coord$Trait  == ">10cm"] <- "Sediment Position"

coord$Trait [coord$Trait  == "Sessile"] <- "Mobility"
coord$Trait [coord$Trait  == "Crawl/creep/climb"] <- "Mobility"
coord$Trait [coord$Trait  == "Burrower"] <- "Mobility"
coord$Trait [coord$Trait  == "Swim"] <- "Mobility"

# Body size
coord$Trait[coord$Trait == "<10"] <- "Maxsize"
coord$Trait[coord$Trait == "500"] <- "Maxsize"
coord$Trait[coord$Trait == "101-200"] <- "Maxsize"
coord$Trait[coord$Trait == "11-20"] <- "Maxsize"
coord$Trait[coord$Trait == "201-500"] <- "Maxsize"
coord$Trait[coord$Trait == "21-100"] <- "Maxsize"

# Longevity
coord$Trait[coord$Trait == "<1"] <- "Longevity"
coord$Trait[coord$Trait == ">10"] <- "Longevity"
coord$Trait[coord$Trait == "1-2"] <- "Longevity"
coord$Trait[coord$Trait == "3-10"] <- "Longevity"

# Larval development
coord$Trait[coord$Trait == "Benthic (Direct)"] <- "Larval Development"
coord$Trait[coord$Trait == "Lecithotrophic"] <- "Larval Development"
coord$Trait[coord$Trait == "Planktotrophic"] <- "Larval Development"

# Feeding mode
coord$Trait[coord$Trait == "Parasite"] <- "Feeding Mode"
coord$Trait[coord$Trait == "Predator"] <- "Feeding Mode"
coord$Trait[coord$Trait == "Scavenger/Opportunist"] <- "Feeding Mode"
coord$Trait[coord$Trait == "Sub-surface Deposit"] <- "Feeding Mode"
coord$Trait[coord$Trait == "Surface Deposit"] <- "Feeding Mode"
coord$Trait[coord$Trait == "Suspension"] <- "Feeding Mode"

# Bioturbation
coord$Trait[coord$Trait == "Diffusive Mixing"] <- "Bioturbation"
coord$Trait[coord$Trait == "Downwards Conveyer"] <- "Bioturbation"
coord$Trait[coord$Trait == "None"] <- "Bioturbation"
coord$Trait[coord$Trait == "Surface Deposition"] <- "Bioturbation"
coord$Trait[coord$Trait == "Upward Conveyor"] <- "Bioturbation"

coord

## Update column order
coord <- coord[,c(7,1:6),]

## Update row order
coord <- coord[match(desiredOrder,coord$Modality), ] 
str(coord)
head(coord)

# Remove NA rows
coord <-coord[complete.cases(coord), ]

## Change colnames for dimensions
names(coord)
colnames(coord)[3] <- "1"
colnames(coord)[4] <- "2"
colnames(coord)[5] <- "3"
colnames(coord)[6] <- "4"
colnames(coord)[7] <- "5"

## Produce table
library(gt)
coord_table <- coord%>%gt()%>%
  fmt_number(columns = vars('1','2','3','4','5'), decimals = 2)%>%
  tab_spanner(
    label = md("**Dimension**"),
    columns = vars('1','2','3','4','5'))%>%
  cols_label(
    Modality = md("**Modality**"),
    Trait = md("**Trait**"),
    '1' = md("**1**"),
    '2' = md("**2**"),
    '3' = md("**3**"),
    '4' = md("**4**"),
    '5' = md("**5**"))#%>%
## Colour cells according to value (https://www.allisonhorst.com/post/2020-03-02-gt-tables-examples/)
#data_color( # Update cell colors...
#  columns = vars('1'), # ...for dose column 
#  colors = scales::col_numeric( # <- bc it's numeric
#    palette = c(
#      "#737373","white","#737373"), # A color scheme (gradient)
#    domain = c(-1,1) # Column scale endpoints
#  )
#)%>%
#  data_color( # Update cell colors...
#    columns = vars('2'), # ...for dose column 
#    colors = scales::col_numeric( # <- bc it's numeric
#      palette = c(
#        "white","#737373"), # A color scheme (gradient)
#      domain = c(-1,1) # Column scale endpoints
#    )
#  )
## View table
coord_table

## Save Response traits table
#CHOOSE
coord_table %>%
  gtsave(
    "OUTPUTS/RESPONSE_TRAITS/rt_pca_coord_table.png")
#OR
coord_table %>%
  gtsave(
    "OUTPUTS/EFFECTS_TRAITS/et_pca_coord_table.png")
#_______________________________________________________________________________
#### RANDOM FOREST: PREPARE DATA ####

## Start with clusteringresults df
names(samclus)

## Get columns in correct order
FaunalCluster=samclus[,c(1,3,2,4)]
#View(FaunalCluster)

## Change names of cols
colnames(FaunalCluster)=c("Sample","lon","lat","cluster")
head(FaunalCluster)

## Number of samples
dim(FaunalCluster)# 18336    4

## Take only the coordinates
FaunalCluster2 <- FaunalCluster[,2:3]
head(FaunalCluster2)
#_______________________________________________________________________________
#### RANDOM FOREST: CREATE RASTER STACK FOR ENV PREDICTOR VARIABLES ####

library(raster)
library(rgdal)

# Hashed out layers not used as covariates (see PrepareRastersxx.R)

bathy <- raster("DATA/Rasters/Final/bathy3.tif")
#chla <- raster("DATA/Rasters/Final/chla3.tif")
cur <- raster("DATA/Rasters/Final/cur3.tif")
gravel <- raster("DATA/Rasters/Mitchell/Predicted_Gravel_Fraction.tif")
#iron <- raster("DATA/Rasters/Final/iron3.tif")
light <- raster("DATA/Rasters/Final/light3.tif")
mud <- raster("DATA/Rasters/Mitchell/Predicted_Mud_Fraction.tif")
#nitr <- raster("DATA/Rasters/Final/nitr3.tif")
oxy <- raster("DATA/Rasters/Final/oxy3.tif")
#phos <- raster("DATA/Rasters/Final/phos3.tif")
phyto <- raster("DATA/Rasters/Final/phyto3.tif")
#prod <- raster("DATA/Rasters/Final/prod3.tif")
sal <- raster("DATA/Rasters/Final/sal3.tif")
#sand <- raster("DATA/Rasters/Mitchell/Predicted_Sand_Fraction.tif")
sil <- raster("DATA/Rasters/Final/sil3.tif")
#spm <- raster("Final/spm3.tif")
spm <- raster("DATA/Rasters/Mitchell/SPM_MEAN.tif")
temp <- raster("DATA/Rasters/Final/temp3.tif")
wov <- raster("DATA/Rasters/Mitchell/Wave_veloc.tif")

# Create raster stack
predictors <- stack(bathy,cur,gravel,light,mud,oxy,phyto,sal,sil,spm,temp,wov)
predictors

# Simple names for predictor variables
names(predictors)=c("Bathymetry","Current","Gravel","Light","Mud","Oxygen","Phyto","Salinity","Silicate","SPM","Temperature","WOV")#"Sand",
names(predictors)

## Plot raster stack
plot(predictors)
#_______________________________________________________________________________
#### RANDOM FOREST: EXTRACT PREDICTOR VARIABLES FROM RASTER STACK ####

head(FaunalCluster)

## Output data for producing confidence layer 
dim(FaunalCluster)
write.csv(FaunalCluster, file = 'C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/OUTPUTS/RESPONSE_TRAITS/FaunalClusterRT.csv')# RESPONSE TRAITS
write.csv(FaunalCluster, file = 'C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/OUTPUTS/EFFECTS_TRAITS/FaunalClusterET.csv')# EFFECTS TRAITS

## Unload extract function from tidyr package (otherwise it won't work)
#.rs.unloadPackage("tidyr")

## Create a df for predictor variables
sdata <- raster::extract(predictors, FaunalCluster[,2:3])
#View(sdata)

##Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(FaunalCluster$Sample,FaunalCluster$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"
head(sdata2)

## Change cols to appropriate type
str(sdata2)
sdata2$Cluster=as.factor(sdata2$Cluster)
#sdata2$Sample=as.character(sdata2$Sample)
#sdata2$AvCur=as.numeric(sdata2$AvCur)
#sdata2$Chla=as.numeric(sdata2$Chla)
#sdata2$Depth=as.numeric(sdata2$Depth)
#$Gravel=as.numeric(sdata2$Gravel)
#sdata2$Mud=as.numeric(sdata2$Mud)
#sdata2$Sal=as.numeric(sdata2$Sal)
#sdata2$Sand=as.numeric(sdata2$Sand)
#sdata2$SPM=as.numeric(sdata2$SPM)
#sdata2$Stress=as.numeric(sdata2$Stress)
#sdata2$Temp=as.numeric(sdata2$Temp)
#sdata2$WOV=as.numeric(sdata2$WOV)

## First check cols of correct type
str(sdata2)
dim(sdata2)# 18336    14
#_______________________________________________________________________________
#### RANDOM FOREST: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Call library
#install.packages("caTools")
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 12 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## Check it was 90% 10% split

## The training set
train = sdata2[ msk,] 
dim(train)#16502    14
head(train)

## Remove station labels for train
train2 =train[,2:14]
#View(train2)

## The test set
test  = sdata2[ !msk,]
dim(test)#1834   14
#View(test)

## Remove station labels for test
test2 =test[,2:14]
#View(test2)
str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata) # 18336    12
dim(train)+dim(test)# 18336    28
#_______________________________________________________________________________
#### RANDOM FOREST: DO MODELLING ####

## Call library
#install.packages("randomForest")
library(randomForest)

## Model
model <- factor(Cluster) ~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV

## Run model
rf2 <- randomForest(model, data=train2,na.action=na.exclude)
#_______________________________________________________________________________
#### RANDOM FOREST: EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####

## Produce plots 
varImpPlot(rf2)

png('OUTPUTS/RESPONSE_TRAITS/variables_affecting_model_response.png') #, height=nrow(pr), width=ncol(pr) RESPONSE TRAITS OR
#png('OUTPUTS/EFFECTS_TRAITS/variables_affecting_model_effects.png') #, height=nrow(pr), width=ncol(pr) EFFECTS TRAITS
varImpPlot(rf2)
dev.off()


preds <- names(rf2$forest$xlevels)

for (i in 1:length(preds)){
  partialPlot(rf2, train2, preds[i], which.class ='1')
  next
}
#_______________________________________________________________________________
#### RANDOM FOREST: EVALUATE THE MODEL PERFORMANCE ####

## Predict cluster group for test set
pred <- predict(rf2, newdata = test2)
table(pred, test2$Cluster)

## We can test the accuracy as follows:
(218+180+198+284+27+60)/ nrow(test2)# RESPONSE TRAITS 52.7%
(198+393+40+59+252+15)/ nrow(test2)# EFFECTS TRAITS 52.1%

## Matches for nearest cluster neighbour(s), based on dendrogram
#(180+218+54+32+198+59+284+49+27+22+60+19)/ nrow(test2)# RESPONSE TRAITS 65.6%

## Confusion matrix plot
#https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot

confusion_matrix <- as.data.frame(table(pred, test2$Cluster))

cm <- ggplot(confusion_matrix, aes(pred,sort(Var2,decreasing = T), fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#737373") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("6","5","4","3","2","1"))

## Save confusion matrix
png('OUTPUTS/RESPONSE_TRAITS/rt_confusion_matrix_plot.png')
#png('OUTPUTS/EFFECTS_TRAITS/et__confusion_matrix_plot.png')
cm
dev.off()
#_______________________________________________________________________________
#### RANDOM FOREST: PRODUCE FULL COVERAGE RASTER ####

##Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)

#_______________________________________________________________________________
#### RANDOM FOREST: OUTPUT RASTER AS TIFF ####

## Save raster
# CHOOSE EITHER:
writeRaster(pr,'DATA/ResponseTraitsCluster.tif',overwrite=TRUE,format = "GTiff")
# OR:
writeRaster(pr,'DATA/EffectsTraitsCluster.tif',overwrite=TRUE,format = "GTiff")
#_______________________________________________________________________________
#### SIMPLE RASTER PLOT (NO BATHY) ####

## Colours in correct order (see results of cluster where nos are matched to groups)
# CHOOSE EITHER
colours <- c('#4575b4', '#00E600', '#91bfdb', '#e0f3f8','#fee090','#fc8d59')# RESPONSE TRAITS OR
#colours <- c('#5ab4ac', '#d8b365', '#c7eae5', '#8c510a','#f6e8c3', '#FF0000')# EFFECTS TRAITS#'#01665e', FF0000

## Plot raster minus legend
plot(pr,col=colours,legend = FALSE)
#plot(pr,col=colours,legend = TRUE)

## Add category legend
legend("bottomright", legend = c("2","1","3","4","5","6"), fill = c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59'))# RESPONSE TRAITS
#legend("bottomright",legend = c("6","3","1","2","4","5"), fill = c('#FF0000', '#c7eae5','#5ab4ac', '#d8b365', '#8c510a', '#f6e8c3'))# EFFECTS TRAITS

#_______________________________________________________________________________
#### 2D RASTER CLUSTER PLOT WITH UNDERLYING BATHY (AS HILLSHADE) ####

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/')

## Load libraries
library(raster)

## Load raster data layers
bathy <- raster("DATA/Rasters/Final/bathy3.tif")
pr = raster('DATA/ResponseTraitsCluster.tif')
#pr = raster('DATA/EffectsTraitsCluster.tif')

## Exaggerate the vertical scale (bathy not very clear unless you do this)
bathy2 <- bathy*bathy

## Generate 'slope' and 'aspect' using the terrain function
slope <- terrain(bathy2,opt='slope')
aspect <- terrain(bathy2,opt='aspect')

## Generate hillshade
hill <- hillShade(slope,aspect,20,0)

# Define a colour palatte (see https://mycolor.space/?hex=%23C2BF5E&sub=1)
# CHOOSE EITHER
colours <- c('#4575b4', '#00E600', '#91bfdb', '#e0f3f8','#fee090','#fc8d59')# RESPONE TRAITS OR
#colours <- c('#5ab4ac', '#d8b365', '#c7eae5', '#8c510a','#f6e8c3', '#FF0000')# EFFECTS TRAITS

## Produce plot
# CHOOSE EITHER
png('OUTPUTS/RESPONSE_TRAITS/raster_response_2Dhillshadebathy.png') # RESPONSE TRAITS OR , width=3000, height=3000, res=400
#png('OUTPUTS/EFFECTS_TRAITS/raster_effects_2Dhillshadebathy.png') # EFFECTS TRAITS
plot(hill,col=grey(1:100/100),legend=FALSE)#,axes=F,box=FALSE
plot(pr, col=colours,add=TRUE,alpha=0.7,axes=F,box=FALSE,legend=FALSE)#no axes
legend(x = 7, y = 53, legend = c("2","1","3","4","5","6"), fill = c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59'),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE
#legend(x = 7, y = 53, legend = c("6","3","1","2","4","5"), fill = c('#FF0000', '#c7eae5','#5ab4ac', '#d8b365', '#8c510a', '#f6e8c3'),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) EFFECTS
dev.off()

#_______________________________________________________________________________
#### PRODUCE SIDE BY SIDE PLOT FOR RESPONSE TRAITS ####

#### 2D RASTER CLUSTER PLOT WITH UNDERLYING BATHY (AS HILLSHADE) ####

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/')

## Load libraries
library(raster)

## Load raster data layers
bathy <- raster("DATA/Rasters/Final/bathy3.tif")
pr = raster('Y:/C8210_OneBenthic/Working_Area/SDM/OUTPUT/ResponseTraitMaxClass.tif') #Anna's version
#pr = raster('DATA/EffectsTraitsCluster.tif')

## Exaggerate the vertical scale (bathy not very clear unless you do this)
bathy2 <- bathy*bathy

## Generate 'slope' and 'aspect' using the terrain function
slope <- terrain(bathy2,opt='slope')
aspect <- terrain(bathy2,opt='aspect')

## Generate hillshade
hill <- hillShade(slope,aspect,20,0)

# Define a colour palatte (see https://mycolor.space/?hex=%23C2BF5E&sub=1)
# CHOOSE EITHER
colours <- c('#4575b4', '#00E600', '#91bfdb', '#e0f3f8','#fee090','#fc8d59')# RESPONE TRAITS OR


## RESPONSE CONFIDENCE PLOT ##
library(RColorBrewer)
prc = raster('Y:/C8210_OneBenthic/Working_Area/SDM/OUTPUT/ResponseTraitConfidence.tif') #Anna's version

cols <- colorBin(palette = brewer.pal(n = 4, name = "Oranges"), bins = c(0,0.10,0.25,0.6), domain = c(0,0.10,0.25,0.6),na.color = "transparent")


plot(hill,col=grey(1:100/100),legend=FALSE)#,axes=F,box=FALSE
plot(pr, col=colours,add=TRUE,alpha=0.7,axes=F,box=FALSE,legend=FALSE)#no axes

#my.palette <- brewer.pal(n = 5, name = "Oranges")
my.palette <- c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603")

cuts=c(0,0.2,0.4,0.6,0.8,1.0) #set breaks,1.0


## SIDE BY SIDE PLOTS WITH a) and b) but on PLOT

png('OUTPUTS/RESPONSE_TRAITS/traits_response_modelplot_confidenceplot.png',width = 30,height = 13.3,units = "cm", res = 600,pointsize = 12)
# 2. Create the plot
line = 1
cex = 1.5
side = 3
adj=-0.15

par(mfrow=1:2)
par(mar= c(2, 3.1, 2, 4)+ 0.2)
plot(pr, col=colours,add=F,alpha=1,axes=T,box=T,legend=F)
legend(x = 8, y = 53.2, legend = c("2","1","3","4","5","6"), fill = c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59'),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE
#
mtext("a)", side=side, line=line, cex=cex, adj=adj)
plot(prc,breaks=cuts,col=my.palette,add=F,alpha=1,axes=T,box=T,legend=F)
legend(x = 4.8, y = 52.5, legend = c("0.0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1.0"), fill = c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE

mtext("b)", side=side, line=line, cex=cex, adj=adj)
dev.off()


#_______________________________________________________________________________
#### PRODUCE SIDE BY SIDE PLOT FOR EFFECTS TRAITS ####

#### 2D RASTER CLUSTER PLOT WITH UNDERLYING BATHY (AS HILLSHADE) ####

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/')

## Load libraries
library(raster)

## Load raster data layers
bathy <- raster("DATA/Rasters/Final/bathy3.tif")
pr = raster('Y:/C8210_OneBenthic/Working_Area/SDM/OUTPUT/EffectsTraitMaxClass.tif') #Anna's version
#pr = raster('DATA/EffectsTraitsCluster.tif')

## Exaggerate the vertical scale (bathy not very clear unless you do this)
bathy2 <- bathy*bathy

## Generate 'slope' and 'aspect' using the terrain function
slope <- terrain(bathy2,opt='slope')
aspect <- terrain(bathy2,opt='aspect')

## Generate hillshade
hill <- hillShade(slope,aspect,20,0)

## Colours
colours <- c('#5ab4ac', '#d8b365', '#c7eae5', '#8c510a','#f6e8c3', '#FF0000')# EFFECTS TRAITS



## EFFECTS CONFIDENCE ##
library(RColorBrewer)
prc = raster('Y:/C8210_OneBenthic/Working_Area/SDM/OUTPUT/EffectsTraitConfidence.tif') #Anna's version

#cols <- colorBin(palette = brewer.pal(n = 4, name = "Oranges"), bins = c(0,0.10,0.25,0.6), domain = c(0,0.10,0.25,0.6),na.color = "transparent")

plot(hill,col=grey(1:100/100),legend=FALSE)#,axes=F,box=FALSE
plot(pr, col=colours,add=TRUE,alpha=0.7,axes=F,box=FALSE,legend=FALSE)#no axes

#my.palette <- brewer.pal(n = 4, name = "Oranges")
my.palette <- c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603")

cuts=c(0,0.2,0.4,0.6,0.8,1.0) #set breaks,1.0

## SIDE BY SIDE PLOTS WITH a) and b) but on PLOT

png('OUTPUTS/EFFECTS_TRAITS/traits_effects_modelplot_confidenceplot.png',width = 30,height = 13.5,units = "cm", res = 600,pointsize = 12)
# 2. Create the plot
line = 1
cex = 1.5
side = 3
adj=-0.15

par(mfrow=1:2)
par(mar= c(2, 3.1, 2, 4)+ 0.2)
plot(pr, col=colours,add=F,alpha=1,axes=T,box=T,legend=F)#
legend(x = 8, y = 53.2, legend = c("6","3","1","2","4","5"), fill = c('#FF0000', '#c7eae5','#5ab4ac', '#d8b365', '#8c510a', '#f6e8c3'),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) EFFECTS
mtext("a)", side=side, line=line, cex=cex, adj=adj)
plot(prc,col=my.palette,add=F,alpha=1,axes=T,box=T,legend=F)
legend(x = 4.8, y = 52.5, legend = c("0.0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1.0"), fill = c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE
mtext("b)", side=side, line=line, cex=cex, adj=adj)

dev.off()
#_______________________________________________________________________________