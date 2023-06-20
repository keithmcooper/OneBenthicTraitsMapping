################################################################################
################################################################################
####                                                                        ####
####      MAPPING MARINE BENTHIC BIOLOGICAL TRAITS TO FACILITATE FUTURE     ####
####      FUTURE SUSTAINABLE DEVELOPMENT (ECOLOGICAL APPLICATIONS)          ####
####                                                                        ####
################################################################################
################################################################################

# This script is split into two parts. Part A deals with the generation of trait 
# cluster groups. Part B deals with the spatial modelling (random forest) of 
# trait cluster groups. Data used in the script is sourced from the OneBenthic
# database using sql queries. For users without access to this database, data
# can be sourced from https://doi.org/10.14466/CefasDataHub.137

################################################################################
####                              PART A                                   #####
####                          Trait Clustering                             #####
################################################################################

## Set working directory
#setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTraitsMapping/')

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
ss.survey_surveyname as surveyname,
s.samplecode,
s.samplelat,
s.samplelong,
wtr.traits_modality,
SUM(wtr.fuzzyscore * sqrt(sqrt(ts.abund))) as traitscore,
SUM  (sqrt(sqrt(ts.abund))) as total_trans_abund,
g.gearcode

FROM samples.sample as s 
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 
INNER JOIN faunal_data.worrmstraits as wtr ON wtr.worrms_aphiaid = w.aphiaid
INNER JOIN gear.gear as g on s.gear_gearcode = g.gearcode
INNER JOIN associations.surveysample as ss on s.samplecode = ss.sample_samplecode
INNER JOIN associations.survey as su on su.surveyname = ss.survey_surveyname

WHERE (wtr.traits_trait = 'morphology               ' OR wtr.traits_trait = 'eggdevelopment           ' OR wtr.traits_trait = 'livinghabit              ' OR wtr.traits_trait = 'sedimentposition         ' OR wtr.traits_trait = 'mobility                 ')
AND (s.gear_gearcode = 'MHN' OR s.gear_gearcode = 'DG' OR s.gear_gearcode = 'VV' OR s.gear_gearcode = 'SM' OR s.gear_gearcode = 'NIOZ' OR s.gear_gearcode = 'BC_0.1' OR s.gear_gearcode = 'C/VV'OR s.gear_gearcode = 'BC')
AND (s.treatment = 'R' or s.treatment IS NULL)
AND s.macrosieve= 1
AND w.include = TRUE
AND w.include_traits = TRUE
AND s.samplelat > 47.92938
AND s.id <= 44362

GROUP BY
ss.survey_surveyname,
s.samplecode,
g.gearcode,
wtr.traits_modality
ORDER by s.samplecode;")

# OR:

##2. SQL select query EFFECTS TRAITS
data = dbGetQuery(pool,
                  "SELECT
ss.survey_surveyname as surveyname,
s.samplecode,
s.samplelat,
s.samplelong,
wtr.traits_modality,
SUM(wtr.fuzzyscore * sqrt(sqrt(ts.abund))) as traitscore,
SUM  (sqrt(sqrt(ts.abund))) as total_trans_abund,
g.gearcode

FROM samples.sample as s 
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 
INNER JOIN faunal_data.worrmstraits as wtr ON wtr.worrms_aphiaid = w.aphiaid
INNER JOIN gear.gear as g on s.gear_gearcode = g.gearcode
INNER JOIN associations.surveysample as ss on s.samplecode = ss.sample_samplecode
INNER JOIN associations.survey as su on su.surveyname = ss.survey_surveyname

WHERE (wtr.traits_trait = 'maxsize                  ' OR wtr.traits_trait = 'longevity                ' OR wtr.traits_trait = 'larvaldevelopment        ' OR wtr.traits_trait = 'feedingmode              ' OR wtr.traits_trait = 'bioturbation             ')
AND (s.gear_gearcode = 'MHN' OR s.gear_gearcode = 'DG' OR s.gear_gearcode = 'VV' OR s.gear_gearcode = 'SM' OR s.gear_gearcode = 'NIOZ' OR s.gear_gearcode = 'BC_0.1' OR s.gear_gearcode = 'C/VV'OR s.gear_gearcode = 'BC')
AND (s.treatment = 'R' or s.treatment IS NULL)
AND s.macrosieve= 1
AND w.include = TRUE
AND w.include_traits = TRUE
AND s.samplelat > 47.92938
AND s.id <= 44362

GROUP BY
ss.survey_surveyname,
s.samplecode,
g.gearcode,
wtr.traits_modality
ORDER by s.samplecode;")

#_______________________________________________________________________________
#### RETRIEVE DATA FROM .csv files (see https://doi.org/10.14466/CefasDataHub.137) #### 

## Set working directory
setwd("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTraitsMapping\\DATA\\FINAL_DOI")

# CHOOSE EITHER:

## 1.Load RESPONSE TRAITS dataset from .csv file
data=read.csv("response_public.csv", header=T,na.strings=c("NA", "-","?","<null>"),
              stringsAsFactors=F,check.names=FALSE)

# OR:

## 2.Load EFFECTS TRAITS dataset from .csv file
data=read.csv("effects_public.csv", header=T,na.strings=c("NA", "-","?","<null>"),
              stringsAsFactors=F,check.names=FALSE)

#_______________________________________________________________________________
## Check data
head(data)

## Check the number of samples
df_uniq <- unique(data$samplecode)
length(df_uniq) #31838

## Create a df for sample codes and gear - this will be used later for modelling
sample_gear <- unique(data[,c(2,8)])
colnames(sample_gear) <- c("sample","gear")

## Calculate traitscore as a percentage of total transformed abundance
data$tperc <- (data$traitscore/data$total_trans_abund)*100

# Drop unwanted columns
data <- data[,c(2,3,4,5,9)]
colnames(data)[5] <- "traitscore"

#_______________________________________________________________________________
#### PREPARE TRAITS DATA FOR CLUSTERING ####

## Examine col names of retrieved data
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
png("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTraitsMapping/OUTPUTS/FIGURE_1.jpg", width = 15, height = 13, units = "cm", res = 800,pointsize = 12) # RESPONSE TRAITS OR
#png("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTraitsMapping/OUTPUTS/FIGURE_3.jpg", width = 15, height = 13, units = "cm", res = 800,pointsize = 12) # EFFECTS TRAITS

## Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(hc)

## Reduce label size of cluster labels 1-5
#par(cex=0.8)
par(cex=0.8)

## Plot dendrogram
plot(hcd, ylab = "Height",leaflab = "none")
#plot(hcd, ylab = "Height")# with labels

## Add circle symbols below each leaf of the dendrogram. Define colours.Useful sites for colours: https://colorbrewer2.org/#type=diverging&scheme=PRGn&n=6 or 
col.circle=c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59') # RESPONE TRAITS https://mycolor.space/?hex=%23BAC25E&sub=1 ##F6906C
#col.circle=c('#FF0000', '#c7eae5','#5ab4ac', '#d8b365', '#8c510a', '#f6e8c3') # EFFECTS TRAITS https://mycolor.space/?hex=%23BAC25E&sub=1

## Add symbols
symbols(1:6, rep(0, 6), circles=rep(1, 6), add=TRUE, inches=.08,fg=col.circle,
        bg=col.circle, xpd=TRUE)#-10

## Add cluster group labels (check with original labels)
axis(1,at=seq(1,6,by=1),labels=c("2","1","3","4","5","6"),pos=-5,cex.axis=0.8,lty = 0,cex.axis=1)# RESPONSE
#axis(1,at=seq(1,6,by=1),labels=c("6","3","1","2","4","5"),pos=-5,cex.axis=0.8,lty = 0,cex.axis=1)# EFFECTS

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
#### PCA (not used in paper) ####
# see https://cran.r-project.org/web/packages/factoextra/readme/README.html.

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
#### PCA CONTRIBUTIONS TABLE (not used in paper)####

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
colnames(FaunalCluster)=c("sample","lon","lat","cluster")
head(FaunalCluster)

## Add in gear type info
head(sample_gear)

## Run this for the response data
resp <- merge(FaunalCluster,sample_gear,by="sample")
resp <- resp[,c(1,2,3,5,4)]## Change order of columns
colnames(resp)[5] <- "resptrait"
head(resp)

## Run this for the effects data
eff <- merge(FaunalCluster,sample_gear,by="sample")
eff <- eff[,c(1,2,3,5,4)]## Change order of columns
colnames(eff)[5] <- "efftrait"
head(eff)

## Merge cluster results for response and effects groups
dataformodelling <- merge(resp,eff,by=c("sample","lon","lat","gear"))
head(dataformodelling)

## save cluster results to csv
write.csv(dataformodelling, file = 'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTraitsMapping\\DATA\\dataformodelling.csv')

## Number of samples
dim(dataformodelling)# 16682

#_______________________________________________________________________________
################################################################################
####                             PART B                                    #####
####                        Random Forest Modelling                        #####
####                      18/06/2023 - Anna Downie                         #####
################################################################################

##### 1. Set up directories, libraries and colour palettes  ####################

## Load required libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(caret)
require(randomForest)
require(pdp)
require(data.table)
require(tmap)
require(sf)
require(raster)
require(ggcorrplot)
require(ggdendro)
require(ggpubr)
require(gridExtra)
require(kableExtra)
require(vtable)
require(caTools)
require(patchwork)


## Set working directory
wkd <- setwd('C:/Users/AD06/Documents/OneBenthic')

## Set auxiliary GIS data file paths (e.g. shoreline)
shoreline <- 'C:/Users/AD06/Documents/Coastlines/Europe/EuropeESRI_high.shp'

## Source algorithm for calculating the Variance Inflation Factor VIF
## The algorithm used in original code is not publicly available but can 
## be replaced here with any algorithm that calculates VIF - save as corvif
source("corvif.R") 


## Colour palette for plots
cpl <- c('#d4ebe7','#cbbcbb','#f5f1f1','#172957','#66afad')
names(cpl) <- c('lt','dbe','lbe','dbl','dt')

## Colour palettes for trait classes
# Colours for Response Traits
classpal.r <- c("#4575b4","#00E600","#91bfdb","#e0f3f8","#fee090", "#fc8d59")
# Colours for Effects Traits
classpal.e <- c("#5ab4ac","#d8b365","#c7eae5","#8c510a","#f6e8c3","#FF0000") 


##### 2. Settings for tables and figures   #####################################

# Table settings for histograms
t1 <- ttheme_minimal(core=list(bg_params = list(fill = '#f5f1f1', col='#172957'),
                               fg_params=list(fontface=1)),
                     colhead=list(fg_params=list(col="#172957", fontface=4L),
                                  bg_params = list(fill = '#f5f1f1', col='#172957')),
                     rowhead=list(fg_params=list(col="#172957", fontface=2L),
                                  bg_params = list(fill = '#f5f1f1', col='#172957')),
                     padding = unit(c(5, 5), "mm"))

# Settings for ggplots
OneB_theme <-
  ggplot2::theme(axis.title.y = element_text(vjust=4,  size=12,colour="black"),
                 axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
                 axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
                 axis.title.x  = element_text(vjust=-4, size=12,colour="black"),
                 strip.background = element_rect(fill=cpl['lbe']),
                 strip.text.x = element_text(size=12, face="bold"),
                 panel.grid.major = element_line(colour=cpl['lbe']),
                 panel.grid.minor = element_line(colour=cpl['lbe']),
                 panel.background = element_rect(fill="white"),
                 plot.margin = ggplot2::margin(0.5, 0.5, 1, 0.5, "cm"))


##### 3. Prepare data   ########################################################


### Read in shoreline for plotting ----
countries <- st_read(shoreline)

### Read in and check response data ----

## Read in point data on traits (including response and effects traits)
## objects is named resptrait as it is the model response variable
resptrait <- as.data.table(read.csv('DATA/ADDITIONAL/dataformodelling.csv'))
rtr_sf <- st_as_sf(resptrait,coords=c('lon','lat'),crs=4326,remove = FALSE)
## Set bounding box for plotting
extbb <- st_bbox(rtr_sf)

## Plot points on map
tm_shape(countries,bbox = extbb) +
  tm_grid(lines = FALSE,n.x = 4, n.y = 3) +
  tm_polygons() +
  tm_shape(rtr_sf) +
  tm_symbols(col="steelblue",size = 0.2)

## List of sampling gears in data
unique(rtr_sf$gear)

# Keeping all 0.1m2 gears (Mini Hamon,Smith-Macintyre, Day and vanveen grabs, 
# nioz box core) and box cores with unknown dimensions, with the expectation 
# they are 0.1m2 

## Combine gears to common types for inclusion as model predictor variable
# Group
rtr_sf <- rtr_sf %>%
  mutate(GearType=case_when(gear %in% c("DG","VV","SM") ~ 'Cup',
                            gear =="MHN" ~'Scoop',
                            gear %in% c("BC_0.1","BC","NIOZ") ~ 'Box Core',
                            gear %in% c("C/VV") ~ 'Core/Cup'))

# Extract and edit columns needed for model response variable data matrix 
# (including both trait categories) 
resp_sf <- rtr_sf %>%
  dplyr::select(c("sample","lon","lat","GearType","resptrait","efftrait")) %>%
  dplyr::rename(SampleCode=sample) %>% # rename sample column to 'SampleCode' 
  mutate_if(is.integer,as.character) # trait classes from integer to character
# Print data summary
summary(resp_sf)

## Plot points by gear type
tm_shape(countries,bbox = extbb) +
  tm_grid(lines = FALSE,n.x = 4, n.y = 3) +
  tm_polygons() +
  tm_shape(resp_sf) +
  tm_symbols(col="GearType",size = 0.2) +
  tm_layout(legend.position=c("right", "bottom"),
            legend.bg.color = 'white',
            legend.bg.alpha = 0.75)


### Read in environmetal rasters ----

## List files
f <- list.files(path='DATA/RASTERS', pattern='.tif$', full.names=T)
## Load rasters into a raster stack
s1 <- stack(f,RAT=F)
## View raster stack to check layers
s1
## Plot raster stack to check layers
# define map extent
extbb <- tmaptools::bb(resp_sf,ext = 1.1)
## Plot layers
tm_shape(s1) + 
  tm_raster(style=rep("kmeans",nlayers(s1)),
            palette="RdYlGn",
            
  ) +
  tm_facets(free.scales = TRUE,ncol = 4) +
  tm_layout(legend.position = c("right", "bottom"),
            panel.show = FALSE,)


# Extract raster data as a dataframe
s1d <- raster::as.data.frame(s1)
s1d <- s1d[complete.cases(s1d),]


### Extract raster values to points ----

## Extract values
env <- raster::extract(s1,resp_sf)
# Check data summary
summary(env)

## Define labels to use for environmental variables
envlab <- c('Depth','Chl-a','Current velocity','Iron','Light','LS-Factor','Nitrate','Oxygen','Phosphate','Phytoplankton','Gravel','Mud','Productivity','Rel. Slope Pos.','Salinity','Silicate','Suspended Matter','Temperature','Wave velocity')
names(envlab) <- colnames(env)
# Add gear type to labels
envlab['GearType'] <- 'Gear type'


### Combine species and environmental data ----

# Names of response columns
rn <- resp_sf %>%
  st_drop_geometry() %>%
  names
# Names of potential predictor columns
re <- env %>%
  colnames()
# Combine and reorder
pshpt <- resp_sf %>%
  cbind(env) %>%
  dplyr::select(rn[-c(2:4)],rn[2:3],all_of(re),rn[4])     

## Traits abd gear as factor
pshpt <- pshpt %>% 
  mutate(across(c(resptrait,efftrait,GearType),factor))

# Keep complete cases
pshpt <- pshpt %>%
  drop_na()
# Check summary
summary(pshpt)


# Predictor variable names including latitude and longitude
lstenv <- ncol(pshpt)-1
fstenv <- lstenv-(length(re)+2)
prnames <- names(pshpt[,fstenv:lstenv])
prnames <- prnames[!prnames=='geometry'] # remove geometry column name

# Define factor variables
facvars <- "GearType"
# Define numeric variables
numvars <- prnames[prnames != facvars]

# Response variable names
rspnames <- c("resptrait","efftrait")

# Define the number of class levels
numclass.r <- nlevels(pshpt$resptrait)
numclass.e <- nlevels(pshpt$efftrait)


### Data exploration ----

## Table of summary statistics
vtable::sumtable(pshpt,simple.kable = TRUE)

## Covariance of environmental variables

# Correlation matrix
corr <- pshpt %>%
  st_drop_geometry() %>%
  dplyr::select(all_of(numvars)) %>%
  cor
# Edit variable names
colnames(corr) <- envlab[colnames(corr)]
colnames(corr)[1:2] <- c('Latitude','Longitude')
rownames(corr) <- envlab[rownames(corr)]
rownames(corr)[1:2] <- c('Latitude','Longitude')
# Plot correlation matrix
corplot <- ggcorrplot(corr, method='circle',type = 'upper',hc.order = TRUE)
corplot


# Correlation visualised as a dedrogram
# Cluster by correlation 
sim.by.hclust <- hclust(dist(corr))
names(sim.by.hclust$labels)[1:2] <- c("samplelat","samplelong")
# Plot dendrogram
cordend <-  ggdendrogram(sim.by.hclust,rotate = TRUE,theme_dendro = TRUE) +
  geom_hline(yintercept = 0.7,linetype=2)
cordend


### Set up objects to save model and validation outputs into ----

# Input response data
RD <- NULL
# An  object for appending all variable importances to
VIs <- NULL
# An object for appending model performance  table data
PSas <- NULL
PScs <- NULL
# Observed vs predicted values across test runs
OvsPs <- NULL
# Partial Response Curve data
PRCs <- NULL

# Tables for collecting all model performance statistics
# Class specific validation statistics
forest.class.res <- data.frame(Name=character(0),Run=character(0),Class=character(0),ClassN=numeric(0),
                               Sens=numeric(0),Spec=numeric(0),BA=numeric(0),
                               stringsAsFactors =F)

# Validation statistics for whole model
forest.class.res.all <- data.frame(Name=character(0),
                                   Run=character(0),
                                   N=numeric(0),
                                   Acc=numeric(0),
                                   NIR=numeric(0),
                                   P=numeric(0),
                                   Kappa=numeric(0),
                                   Q=numeric(0),
                                   A=numeric(0),
                                   stringsAsFactors =F)

# Combined validations statistic
forest.res.mat.all <- data.frame(Name=character(0),ModRun = numeric(),Comb=character(),
                                 Pred=character(),Obs=character(),Vals=numeric(),
                                 stringsAsFactors = FALSE)




##### 4. Response traits model   ########################################################

###  Data prep ----

## Define response variable as response traits
tax = 'resptrait'

## Select all environmental variables to include as potential predictors and
## save into a temporary dataset
cols <- c(tax,prnames)
nrow(pshpt)
sdata <- pshpt %>%
  st_drop_geometry() %>%
  dplyr::select(all_of(cols)) %>%
  dplyr::rename(resp=1) %>% # Rename the response column
  drop_na() # Include only complete cases
# Data summary
summary(sdata)
#levels(sdata$resp)[which.min(summary(sdata$resp))]

## Add data to model data collator
RD <- RD %>%
  rbind(data.table(Tax=tax,Metric='Class',value=sdata[,'resp']))


## Histogram of input data
mxy <- max(summary(sdata$resp)) # Find max for y
indat <- ggplot(sdata,(aes(x=resp))) +
  geom_bar(fill='#66afad', col='#172957') +
  geom_text(inherit.aes=FALSE,aes(label = paste('N =',nrow(sdata)), x = 1,
                                  y = mxy),position = position_nudge(y=mxy/10),size=5,hjust=0) +
  xlab('Class') +
  ylab('Number of samples') +
  theme(axis.title.y = element_text(vjust=0.5,size=12,colour="black"),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_text(vjust=-1, size=12,colour="black"),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="#f5f1f1"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"),
        plot.margin = ggplot2::margin(3,3,8,3))
indat


## Plot model data on a map
rtr_sf <- st_as_sf(sdata[,1:3],coords=c('lon','lat'),crs=CRS('+proj=longlat +datum=WGS84 +no_defs'))

asg.ab <- tm_shape(countries,bbox = extbb) +
  tm_grid(lines = FALSE,n.x = 4, n.y = 3) +
  tm_layout(title = "Response trait group" ,
            legend.position = c("right", "bottom"),
            legend.title.fontface = 3,
            legend.bg.color = 'white',
            legend.width = 0.5) +
  tm_polygons() +
  tm_shape(rtr_sf) +
  tm_bubbles(col='resp',
             palette = classpal.r,
             size=0.2,
             alpha=0.5,
             border.col = "black", border.alpha = .5, 
             scale = 1.1,
             title.col = "Group")

asg.ab

### Preliminary full model with all potential predictors ----

## Build model with all variables
prelRF <- randomForest(resp~.,
                       data=sdata, 
                       ntrees=1000,
                       replace=FALSE,
                       importance=TRUE,
                       nPerm=5)
prelRF

## Extract variabe importance
full.importance <- data.table(Predictor=rownames(prelRF$importance),prelRF$importance)
full.importance <- full.importance[order(full.importance[,MeanDecreaseGini],decreasing=T),]
full.importance


## Plot Partial Dependence Curve

# Objects for storing data
plotdata <- NULL
plotdata2 <- NULL
predselnf <- numvars # define numeric variables

# Loop through response levels and predictor variables to produce and save 
# partial dependence curve data for each case
require(pdp)

# Extract plot values for continuous variables
for (j in 1:length(predselnf)) {
  
  for (c in levels(sdata$resp)) {
    
    pdata <- partial(prelRF,pred.var = predselnf[j],which.class = c,
                     plot = FALSE,train=sdata,grid.resolution=100,prob = TRUE)
    predname <- predselnf[j]
    print(predname)
    temp <- data.frame(predvar=predselnf[j],class=c,x=pdata[[1]],y=pdata[[2]])
    plotdata <- rbind(plotdata,temp)
    
  }
  
}

# Round values
plotdata.r <- plotdata %>%
  mutate(across(y, round, 1))

# Extract plot values for categorical variables
for (j in 1:length(facvars)){
  
  for (c in levels(sdata$resp)) {
    
    pdata <- partial(prelRF,pred.var = facvars[j],which.class = c,
                     plot = FALSE,train=sdata,grid.resolution=100,prob = TRUE)
    predname <- facvars[j]
    temp <- data.frame(predvar=facvars[j],class=c,x=pdata[[1]],y=pdata[[2]])
    plotdata2 <- rbind(plotdata2,temp)
    
  }
  
}
# Round values
plotdata2 <- plotdata2 %>%
  mutate(across(where(is.numeric), round, 1))


## Create list object to store individual plots
fullRP.list <- list()

## Populate list object with plots for continuous variables
for (i in predselnf) {
  
  fullRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y,col=class)) +
    geom_smooth(size=0.8,se=FALSE,span = 0.3,show.legend = FALSE) +
    facet_wrap(~ predvar,scales = "free_x", ncol=3,labeller = labeller(predvar=envlab)) +
    scale_colour_manual(values = classpal.r) +
    ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
    theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
          axis.title.x  = element_blank(),
          plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
          strip.background = element_rect(fill="grey90"),
          strip.text.x = element_text(size=12, face="bold"),
          panel.grid.major = element_line(colour="grey80"),
          panel.grid.minor = element_line(colour="grey80"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=12,colour="black"))
}

# Add plot for gear type
fullRP.list[['GearType']] <-  ggplot(plotdata2, aes(x=x,y=y,fill=class)) +
  geom_bar(stat = 'identity',position = "dodge") +
  scale_fill_manual(values = classpal.r) +
  facet_wrap(~ predvar,scales = "free_x", ncol=3,labeller = labeller(predvar=envlab)) +
  ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
  theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_blank(),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey90"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"))

# Create plot layout
fullRP <- grid.arrange(grobs=fullRP.list,ncol=4)
# View plot
fullRP


## Select best performing uncorrelated variables to keep for final model

# Correlation matrix with variables in order of full model importance
# Variable list
vl <- full.importance[[1]]
# Numeriv variables only
vl <- vl[vl %in% numvars]
cr <- sdata %>%
  dplyr::select(all_of(vl)) %>%
  cor()

# Remove variables correlated to a higher importance variable
for(j in 1:length(cr[1,])){
  if (j == 1){
    pl <- c(names(cr[j,][1]),names( cr[j,][sqrt((cr[j,])^2)<0.8]))
    pl1 <- pl
  } else if (names(cr[j,])[j] %in% pl1){
    rem <- names(cr[j,-c(1:j)][sqrt((cr[j,-c(1:j)])^2)>0.8])
    if (length(rem) != 0L){  
      pl <- pl[!pl %in% rem]
    }
  }
  next
}

# Table of kept variables and their Variance Inflation Factors (vif)
crval <- as.data.frame(pl)
crval[2] <- sdata %>%
  dplyr::select(all_of(pl)) %>%
  corvif()
crval


## Same Without Coordinates
# Variable list
vl2 <- vl[!vl %in% c("lon", "lat")]
# Correlation matrix
cr <- sdata %>%
  dplyr::select(all_of(vl2)) %>%
  cor()
# Remove variables correlated to a higher importance variable
for(j in 1:length(cr[1,])){
  if (j == 1){
    pl <- c(names(cr[j,][1]),names( cr[j,][sqrt((cr[j,])^2)<0.8]))
    pl1 <- pl
  } else if (names(cr[j,])[j] %in% pl1){
    rem <- names(cr[j,-c(1:j)][sqrt((cr[j,-c(1:j)])^2)>0.8])
    if (length(rem) != 0L){  
      pl <- pl[!pl %in% rem]
    }
  }
  next
}

# Table of kept variables and their Variance Inflation Factors (vif)
crval2 <- as.data.frame(pl)
crval2[2] <- sdata %>%
  dplyr::select(all_of(pl)) %>%
  corvif()
crval2

### Choose the set of variables to use in model ----

## Set with coordinates
predsel <- crval[[1]]

# List of all variables including response
clms <- c(names(sdata)[1],predsel,facvars)

# Data to use in model
mdata <- sdata %>%
  dplyr::select(all_of(clms))
summary(mdata)

## Full data random forest with selected variables
selRF <- randomForest(resp~.,
                      data=mdata, 
                      ntrees=1000,
                      replace=FALSE,
                      importance=TRUE,
                      nPerm=5)
selRF


## Choosing set without coordinates
predsel <- crval2[[1]]

# List of all variables including response
clms <- c(names(sdata)[1],predsel,facvars)

# Data to use in model
mdata <- sdata %>%
  dplyr::select(all_of(clms))
summary(mdata)

## Full data random forest with selected variables
selRF <- randomForest(resp~.,
                      data=mdata, 
                      ntrees=1000,
                      replace=FALSE,
                      importance=TRUE,
                      nPerm=5)
selRF

## Not much difference in the two models so keeping the one without coordinates

sel.importance <- data.table(Predictor=rownames(selRF$importance),selRF$importance)
sel.importance <- sel.importance[order(sel.importance[,2],decreasing=T),]
sel.importance

## Plot response curves to confirm variables are sensible
# Plot Partial Dependence
plotdata <- NULL
plotdata2 <- NULL

require(pdp)

# Extract plot values
for (j in 1:length(predsel)) {
  
  for (c in levels(sdata$resp)) {
    
    pdata <- partial(selRF,pred.var = predsel[j],which.class = c,
                     plot = FALSE,train=sdata,grid.resolution=100,prob = TRUE)
    predname <- predsel[j]
    print(predname)
    temp <- data.frame(predvar=predsel[j],class=c,x=pdata[[1]],y=pdata[[2]])
    plotdata <- rbind(plotdata,temp)
    
  }
  
}

plotdata.r <- plotdata %>%
  mutate(across(y, round, 1))

for (j in 1:length(facvars)){
  
  for (c in levels(sdata$resp)) {
    
    pdata <- partial(selRF,pred.var = facvars[j],which.class = c,
                     plot = FALSE,train=sdata,grid.resolution=100,prob = TRUE)
    predname <- facvars[j]
    temp <- data.frame(predvar=facvars[j],class=c,x=pdata[[1]],y=pdata[[2]])
    plotdata2 <- rbind(plotdata2,temp)
    
  }
  
}

plotdata2 <- plotdata2 %>%
  mutate(across(where(is.numeric), round, 1))

# Produce plots in a list to combine later
selRP.list <- list()

for (i in predsel) {
  
  selRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y,col=class)) +
    geom_smooth(size=0.8,se=FALSE,span = 0.3,show.legend = FALSE) +
    facet_wrap(~ predvar,scales = "free_x", ncol=3,labeller = labeller(predvar=envlab)) +
    scale_colour_manual(values = classpal.r) +
    ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
    theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
          axis.title.x  = element_blank(),
          plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
          strip.background = element_rect(fill="grey90"),
          strip.text.x = element_text(size=12, face="bold"),
          panel.grid.major = element_line(colour="grey80"),
          panel.grid.minor = element_line(colour="grey80"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=12,colour="black"))
}

selRP.list[['GearType']] <-  ggplot(plotdata2, aes(x=x,y=y,fill=class)) +
  geom_bar(stat = 'identity',position = "dodge",show.legend = FALSE) +
  scale_fill_manual(values = classpal.r) +
  facet_wrap(~ predvar,scales = "free_x", ncol=3,labeller = labeller(predvar=envlab)) +
  ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
  theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black",angle = 0),
        axis.title.x  = element_blank(),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey90"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"))

selRP.list[['Legend']] <- as_ggplot(get_legend(ggplot(plotdata[plotdata$predvar== "bathy3",],aes(x=x,y=y,col=class)) +
                                                 geom_smooth(size=0.8,se=FALSE,span = 0.3) +
                                                 scale_colour_manual(values = classpal.r,name='Response trait class',
                                                                     guide=guide_legend(title.position = 'top')) +
                                                 theme(legend.direction = 'horizontal',
                                                       legend.key = element_blank())))


selRP <- grid.arrange(grobs=selRP.list,ncol=4)



#### After checking plots continuing with the variable set without coordinates #########

###  Cross-validation ----

## Insert number of cross-validation runs required
nruns <- 10 

## Determine names of predictors to be used 
preds <-   c(crval2[[1]],facvars)
predselnf <- predsel[!predsel %in% facvars]

# Set up lists for looping through
train.sets <- list()
test.sets <- list()

# Split for 10 random subsets (list of row numbers)
trainIndex <- createDataPartition(mdata$resp, p = .75,
                                  times = nruns)
# Create train and tests sets
for (j in 1:10){
  
  train.sets[[j]] <- mdata[trainIndex[[j]],]
  test.sets[[j]] <- mdata[-trainIndex[[j]],]
  
  next}


## Run model

# List objects to store outputs
ffs <- list()
imps <- list()
res <- list()
# Tables for collecting all model performance statistics
# Class specific validation statistics
class.res <- data.frame(Name=character(0),Run=character(0),Class=character(0),ClassN=character(0),
                        Sens=numeric(0),Spec=numeric(0),BA=numeric(0),
                        stringsAsFactors =F)

# Validation statistics for whole model
class.res.all <- data.frame(Name=character(0),
                            Run=character(0),
                            N=character(0),
                            Acc=numeric(0),
                            NIR=numeric(0),
                            P=numeric(0),
                            Kappa=numeric(0),
                            Q=numeric(0),
                            A=numeric(0),
                            stringsAsFactors =F)

# Combined validations statistic
res.mat.all <- data.frame(Name=character(0),ModRun = numeric(),Comb=character(),
                          Pred=character(),Obs=character(),Vals=numeric(),
                          stringsAsFactors = FALSE)

plotdata <- NULL
plotdata2 <- NULL

for (j in 1:10){
  
  train <- train.sets[[j]]
  test <- test.sets[[j]]
  
  ffs[[j]] <- randomForest(resp ~.,data=train,
                           ntree=500, replace=FALSE,importance=T, keep.forest= T)
  
  
  results <- as.data.frame(rownames(test))
  results$actual <- test[[1]]
  # Predict class with model i
  results$predicted <- as.data.table(predict(ffs[j],test))[,V1]
  names(results) <- c("id", "actual", "predicted")
  
  # Calculate confusion matrix for predictions by model i
  results.matrix <- confusionMatrix(results$predicted, results$actual)
  results.matrix
  
  # Get the number of objects predicted into each class by model i
  temp.pred <- NULL
  s1dnrow <- nrow(s1d)
  rs <- data.frame(start=c(1,ceiling(s1dnrow/3),2*ceiling(s1dnrow/3)),
                   end=c(ceiling(s1dnrow/3)-1,2*ceiling(s1dnrow/3)-1,s1dnrow))
  for (e in 1:nrow(rs)){
    
    ts1d <- s1d[rs[e,'start']:rs[e,'end'],]
    ts1d$GearType <- factor('Scoop',levels = levels(mdata$GearType))
    temp.pred <- rbind(temp.pred,as.data.table(predict(ffs[j],ts1d,'response')))
    
  } 
  
  pctObj <- as.data.frame(summary(temp.pred[,V1])/sum(summary(temp.pred[,V1])))
  
  # Get values from confusion matrix for model i into a data frame for future plotting
  rnx <- 1
  res.mat.0 <- data.frame(Name=character(), ModRun = numeric(),Comb=character(),
                          Pred=character(),Obs=character(),Vals=numeric(),
                          ValsP=numeric(),
                          stringsAsFactors = FALSE)
  
  for (r in 1:dim(results.matrix$table)[1]){
    
    for (o in 1:dim(results.matrix$table)[2]){
      rn <- rnx+o-1
      res.mat.0[rn,2] <- j
      labt <- paste(r,o,sep="") 
      res.mat.0[rn,3] <- labt
      res.mat.0[rn,4] <- levels(train[[1]])[r]
      res.mat.0[rn,5] <- levels(train[[1]])[o]
      res.mat.0[rn,6] <- results.matrix$table[r,o]
      next
    }
    
    CurSel <- res.mat.0[(rn-o+1):rn,]
    PrVals <- res.mat.0[(rn-o+1):rn,6]
    totPrCL <- sum(res.mat.0[(rn-o+1):rn,6])
    
    res.mat.0[(rn-o+1):rn,7] <- round((PrVals/totPrCL) * pctObj [[1]][r],3)
    CurSel <- res.mat.0[(rn-o+1):rn,]
    
    labt <- paste(r,'P',sep="") 
    rn <- rn+1
    res.mat.0[rn,2] <- j
    res.mat.0[rn,3] <- labt
    res.mat.0[rn,4] <- levels(train[[1]])[r]
    res.mat.0[rn,5] <- 'P'
    res.mat.0[rn,6] <- round(100*(results.matrix$table[r,r]/rowSums(results.matrix$table)[r]),1)
    res.mat.0[rn,7] <- round(100*(CurSel[CurSel[4]==CurSel[5],7]/sum(CurSel[7])),3)
    
    rnx <- rnx+dim(results.matrix$table)[2]+1
    next
  }
  
  
  
  for (o in 1:dim(results.matrix$table)[2]){
    ObsClassVal <- res.mat.0[res.mat.0[2]==j & res.mat.0[5]==levels(train[[1]])[o] ,]
    rn <- rnx+o-1
    res.mat.0[rn,2] <- j
    labt <- paste('U',o,sep="") 
    res.mat.0[rn,3] <- labt
    res.mat.0[rn,4] <- 'U'
    res.mat.0[rn,5] <- levels(train[[1]])[o]
    res.mat.0[rn,6] <- round(100*(results.matrix$table[o,o]/colSums(results.matrix$table)[o]),1)
    res.mat.0[rn,7] <- round(100*(ObsClassVal[ObsClassVal[4]==ObsClassVal[5],7]/sum(ObsClassVal[7])),1)
    next
  }
  
  RunValAll <- res.mat.0[res.mat.0[2]==j,]
  RunValAll$Pred <- factor(RunValAll$Pred,levels = c(levels(train[[1]]),"U","P"))
  RunValAll$Obs <- factor(RunValAll$Obs,levels = c(levels(train[[1]]),"U","P"))
  RunValAll$ValsP[is.nan(RunValAll$ValsP) & RunValAll$Pred !="U" & RunValAll$Pred !="P"] <- 0
  RunValNum <- RunValAll[RunValAll[4]!="U" & RunValAll[5]!="P",]
  RunValCor <- RunValNum[RunValNum[4]==RunValNum[5],]
  Psum <- aggregate(ValsP~Pred,data=RunValNum,sum,na.action=na.pass)
  Osum <- aggregate(ValsP~Obs,data=RunValNum,sum,na.action=na.pass)
  
  res.mat.0[rn+1,2] <- j
  labt <- paste('U','P',sep="") 
  res.mat.0[rn+1,3] <- labt
  res.mat.0[rn+1,4] <- 'U'
  res.mat.0[rn+1,5] <- 'P'
  res.mat.0[rn+1,6] <- round(100*results.matrix[[3]][[1]],1)
  res.mat.0[rn+1,7] <- round(100*(sum(RunValCor[7])/sum(RunValNum[7])),3)
  
  res.mat.all <- rbind(res.mat.all,res.mat.0)
  
  require(matrixStats)
  sum(2*rowMins(as.matrix(cbind(Osum[2]-RunValCor[7],Psum[2]-RunValCor[7]))))/2
  
  # Get overall accuracy measures for model validation run i
  class.res.all[j,2] <- j
  class.res.all[j,3] <- nrow(test)
  class.res.all[j,4] <- results.matrix[[3]][[1]]
  class.res.all[j,5] <- results.matrix[[3]][[5]]
  class.res.all[j,6] <- results.matrix[[3]][[6]]
  class.res.all[j,7] <- results.matrix[[3]][[2]]
  class.res.all[j,8] <- sum(abs(Osum[2]-Psum[2]))/2
  class.res.all[j,9] <- sum(2*rowMins(as.matrix(cbind(Osum[2]-RunValCor[7],
                                                      Psum[2]-RunValCor[7]))))/2
  
  class.res.0 <- data.frame(Name=character(0),Run=character(0),Class=character(0),ClassN=character(0),
                            Sens=numeric(0),Spec=numeric(0),BA=numeric(0),
                            stringsAsFactors =F)
  
  # Get class-specific accuracy measures for model validation run i
  for (i in 1:numclass.r){
    
    class.res.0[i,2] <- j
    class.res.0[i,3] <- row.names(as.data.frame(results.matrix[[4]]))[i]
    class.res.0[i,4] <- sum(results.matrix$table[,i])
    class.res.0[i,5] <- as.data.frame(results.matrix[[4]])[i,1]
    class.res.0[i,6] <- as.data.frame(results.matrix[[4]])[i,2]
    class.res.0[i,7] <- as.data.frame(results.matrix[[4]])[i,11]
    
    next
  }
  
  class.res <- rbind(class.res,class.res.0)
  
  
  ## Store Validation Results in main tables
  class.res$Name <- tax
  class.res
  
  
  class.res.all$Name <- tax
  class.res.all
  
  
  res.mat.all$Name <- tax
  res.mat.all
  
  
  imps[[j]] <- list(round(randomForest::importance(ffs[[j]]), 2))
  
  require(pdp)
  
  for (p in 1:length(predselnf)) {
    
    for (c in levels(mdata$resp)) {
      
      pdata <- partial(ffs[[j]],pred.var = predselnf[p],which.class = c,
                       plot = FALSE,train=mdata,grid.resolution=50,prob = TRUE)
      predname <- predselnf[p]
      temp <- data.frame(Name='resptrait',predvar=predselnf[p],class=c,x=pdata[[1]],y=pdata[[2]])
      plotdata <- rbind(plotdata,temp)
      
    }
    
  }
  
  if (!is.null(facvars)) {
    
    for (f in 1:length(facvars)){
      
      for (c in levels(mdata$resp)) {
        
        pdata <- partial(ffs[[j]],pred.var = facvars[f],which.class = c,
                         plot = FALSE,train=mdata,grid.resolution=50,prob = TRUE)
        predname <- facvars[j]
        temp <- data.frame(Name='resptrait',predvar=facvars[f],class=c,x=pdata[[1]],y=pdata[[2]])
        plotdata2 <- rbind(plotdata2,temp)
        
      }
      
    }
    
    
  }
  
  
  next
}

## Add cross-validation results to tables collecting all values
forest.class.res <- rbind(forest.class.res,class.res)
forest.class.res.all <- rbind(forest.class.res.all,class.res.all)  
forest.res.mat.all <- rbind(forest.res.mat.all,res.mat.all)


### Add confusion matrix to collector list ----
CMTX <-  as.data.table(res.mat.all)

## Create and add performance tables to list ----

## Separate the main matrix and producers, 
## users and overall accuracy into different data frames
res.mat.num <- res.mat.all[res.mat.all[4]!="U" & res.mat.all[5]!="P",]
res.mat.num <- transform(res.mat.num, 
                         Pred = factor(Pred, levels=levels(train[[1]])), 
                         
                         Obs = factor(Obs, levels=levels(train[[1]])))
res.mat.up <- res.mat.all[res.mat.all[4]=="U" | res.mat.all[5]=="P",]
res.mat.up <- transform(res.mat.up, 
                        Pred = factor(Pred, levels = c(levels(train[[1]]),"U","P")),
                        Obs = factor(Obs, levels = c(levels(train[[1]]),"U","P")))
res.mat.up[res.mat.up$Comb=="UP","Comb"] <- "OA"
highlights <- res.mat.num[res.mat.num[4]==res.mat.num[5],][1:numclass.r,4:6]


## Compare between classes
# Rename classes for class specific results for consistency
str(class.res)
class.res$Class <- as.factor(class.res$Class)
levels(class.res[[3]])
levels(class.res[[3]]) <- levels(train[[1]])
class.res$Run <- as.numeric(class.res$Run)
class.res$ClassN <- as.numeric(class.res$ClassN)

# Calculate overall true positives and sensitivity for each model run
TP <- aggregate(Vals~ModRun,data=res.mat.num[res.mat.num$Pred == res.mat.num$Obs,],sum)
SEN <- TP[2]/aggregate(Vals~ModRun,data=res.mat.num,sum)[2]

# Calculate overall true negatives and specififity for each model run
negs <- rep(0,length(unique(res.mat.num[1])[[1]]))
tnegs <- rep(0,length(unique(res.mat.num[1])[[1]]))

for (i in unique(res.mat.num[4])[[1]]){
  
  inegs <- aggregate(Vals~ModRun,data=res.mat.num[res.mat.num$Obs != i,],sum)
  negs <- negs+inegs[[2]]
  itnegs <- aggregate(Vals~ModRun,data=res.mat.num[res.mat.num$Obs != i & res.mat.num$Pred != i ,],sum)
  tnegs <- tnegs+itnegs[[2]]
  
  
  next
}

SPE <- tnegs/negs

# Calculate overall balanced accuracy for each model run   
BA <- (SEN+SPE)/2

# Combine values to a matrix
classrestot <- data.frame('resptrait',1:length(SPE),"Overall",nrow(train),SEN,SPE,BA)
names(classrestot) <- names(class.res)

# Add overall accuracy values to class specific table
classvalB <- rbind(class.res,classrestot)

classvalB
str(classvalB)

# Calculate averages and standard deviations for validation statistics across model runs
cavevalsB <- aggregate(x = classvalB[4:7], by = list(classvalB$Class), FUN = "mean")
csdvalsB <- aggregate(x = classvalB[4:7], by = list(classvalB$Class), FUN = "sd")

# Combine values in a table
BCvalsT <- data.frame(Name=cavevalsB[1],
                      N=cavevalsB[2],
                      SENSmean=round(cavevalsB[3],2),
                      SENSsd=round(csdvalsB[3],2),
                      SPECmean=round(cavevalsB[4],2),
                      SPECsd=round(csdvalsB[4],2),
                      BAmean=round(cavevalsB[5],2),
                      BAsd=round(csdvalsB[5],2))
# Rename columns
names(BCvalsT) <- c("Name","N","SENSmean","SENSsd","SPECmean","SPECsd",
                    "BAmean","BAsd")
# Print table
BCvalsT

asg.cl.perf <- data.table(Cluster=BCvalsT$Name,
                          N =BCvalsT$N,
                          'Sensitivity'= paste(BCvalsT$SENSmean, '\u00B1',BCvalsT$SENSsd),
                          'Specificity' =  paste(BCvalsT$SPECmean, '\u00B1',BCvalsT$SPECsd),
                          'Balanced Accuracy'= paste(BCvalsT$BAmean, '\u00B1',BCvalsT$BAsd))



## Add to the performance data table collator list
PScs <- rbind(PScs,data.table(Tax=tax,asg.cl.perf[-7, ]))

### Plot the class specific accuracies ----

str(classvalB)

require(reshape2)
plotterBC <-melt(classvalB,id.vars=c(3),measure.vars=c(5:7))
str(plotterBC)
levels(plotterBC$variable)
levels(plotterBC$variable) <- c("Sensitivity","Specificity","Balanced Accuracy")

require(ggplot2)
ssbp <- ggplot(plotterBC,(aes(x=Class,y=value,fill=Class))) +
  geom_boxplot(width=0.7,position=position_dodge(width=0.71)) +
  scale_fill_manual(values=c(classpal.r,'#000000')) +
  facet_wrap(~ variable) +
  ylim(c(0,1)) +
  ggtitle("Accuracy Statistics Per Class") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(angle=90,vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_blank(),
        plot.title =  element_text(size=14,colour="black", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey40"),
        strip.text.x = element_text(size=14, face="bold",colour="white"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"))


### Look at the Whole classification ----
class.res.all
classvalB.all <- data.frame(class.res.all,Sens=SEN[[1]],Spec=SPE,BA=BA[[1]])
classvalB.all$N <- as.numeric(classvalB.all$N)

# Calculate averages and standard deviations for validation statistics
callavevalsB <- colMeans(classvalB.all[,4:12])
callsdvalsB <- colSds(as.matrix(classvalB.all[,4:12]))

# Combine values in a table
BCallvalsT <- data.frame(Accmean=round(callavevalsB[1],2),
                         Accsd=round(callsdvalsB[1],2),
                         Pmean=round(callavevalsB[3],2),
                         Psd=round(callsdvalsB[3],2),
                         Kmean=round(callavevalsB[4],2),
                         Ksd=round(callsdvalsB[4],2),
                         Qmean=round(callavevalsB[5],2),
                         Qsd=round(callsdvalsB[5],2),
                         Amean=round(callavevalsB[6],2),
                         Asd=round(callsdvalsB[6],2),
                         Sensmean=round(callavevalsB[7],2),
                         Senssd=round(callsdvalsB[7],2),
                         Specmean=round(callavevalsB[8],2),
                         Specsd=round(callsdvalsB[8],2),
                         BAmean=round(callavevalsB[9],2),
                         BAsd=round(callsdvalsB[9],2))

# Rename columns
names(BCallvalsT) <- c("Accmean","Accsd","Pmean","Psd","Kmean","Ksd",
                       "Qmean","Qsd","Amean","Asd","Sensmean","Senssd",
                       "Specmean","Specsd","BAmean","BAsd")


# Print table
BCallvalsT

asg.perf <- data.table(N = nrow(train),
                       'Sensitivity'= paste(BCallvalsT$Sensmean, '\u00B1',BCallvalsT$Senssd),
                       'Specificity' =  paste(BCallvalsT$Specmean, '\u00B1',BCallvalsT$Specsd),
                       'Kappa' = paste(BCallvalsT$Kmean, '\u00B1',BCallvalsT$Ksd) ,
                       'Balanced Accuracy'= paste(BCallvalsT$BAmean, '\u00B1',BCallvalsT$BAsd),
                       'Quantity Disagreement'=paste(BCallvalsT$Qmean, '\u00B1',BCallvalsT$Qsd),
                       'Allocation Disagreement'=paste(BCallvalsT$Amean, '\u00B1',BCallvalsT$Asd))


## Add to Performance stats collator list
PSas <- rbind(PSas,data.table(Tax=tax,asg.perf[, data.table(t(.SD), keep.rownames=TRUE),]))


#### Observed vs. predicted plot ----

## Plot of confusion matrix and accuracies

## Define the plot of main confusion matrix
## Main plot
p <-  ggplot(res.mat.num,(aes(x=Obs,y=Vals))) +
  geom_boxplot(width=0.7,fill=cpl['dt']) +
  facet_grid(Pred~Obs, scales = "free") +
  ylab("PREDICTED\n ") +
  ggtitle("OBSERVED\n ") + 
  OneB_theme +
  theme(axis.title.y = element_text(size=12,colour="black", face = "bold",vjust=2),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_blank(),
        axis.title.x  = element_blank(),
        plot.title =  element_text(size=12,colour="black", face = "bold",vjust=0.7),
        strip.text.x = element_text(size=12, face="bold",colour="white"),
        strip.text.y = element_text(size=12, face="bold",colour="white"),
        strip.background = element_rect(fill="grey30"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",color="grey30", 
                                        linewidth =0.5, linetype="solid"),
        #plot.margin = unit(c(0.5, 0, 0,0.7), "cm")
  )

p

## Add shading to diagonal
fillcol <- "grey80" 
p1 <- p + geom_rect(data=highlights,aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), 
                    fill=fillcol, alpha=0.3) +
  geom_boxplot(data=res.mat.num,aes(x=Obs,y=Vals),width=0.7,fill=cpl['dt'])
p1


### Define plot of the ranges of users accuracies
p2 <-  ggplot(res.mat.up[grep("U",res.mat.up$Comb),],(aes(x=Comb,y=Vals))) +
  geom_boxplot(width=0.7,position=position_dodge(width=0.71),fill=cpl['dt']) +
  expand_limits(y=c(0,100)) +
  facet_wrap(~Comb, nrow = 1,scales = "free_x") +
  ylab("USERS \n ACCURACY %") +
  ggtitle("  ") + 
  theme(axis.title.y = element_text(size=12,colour="black", face = "bold",vjust=1.7),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title =  element_text(size=18,colour="black", face = "bold",vjust=2),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",color="grey30", size=0.5, linetype="solid"),
        #plot.margin = unit(c(0.5, 0.8, 0, 0.5), "cm")
  )
p2

### Define plot of the ranges of producers accuracies
p3 <-  ggplot(res.mat.up[grep("U",res.mat.up$Comb),],(aes(x=Comb,y=Vals))) +
  geom_boxplot(width=0.7,position=position_dodge(width=0.71),fill=cpl['dt']) +
  expand_limits(y=c(0,100)) +
  facet_wrap(~Comb, ncol=1,scales = "free_x") +
  ggtitle("PRODUCER'S \n ACCURACY %") + 
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title =  element_text(size=12,colour="black", face = "bold",vjust=0),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",color="grey30", 
                                        size=0.5, linetype="solid"),
        #panel.spacing = unit(c(0, 0, 0, 0), "cm"),
        #plot.margin = unit(c(0.4, 0.3, 0.1, 0.5), "cm")
  )
p3

### Define plot of the ranges of overall accuracy
p4 <-  ggplot(res.mat.up[grep("OA",res.mat.up$Comb),],(aes(x=Comb,y=Vals))) +
  geom_boxplot(width=0.7,position=position_dodge(width=0.71),fill=cpl['dt']) +
  expand_limits(y=c(0,100)) +
  facet_wrap(~Comb, ncol=1,scales = "free_x") +
  ggtitle("OVERALL \n ACCURACY %") + 
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title =  element_text(size=12,colour="black", face = "bold",vjust=0),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",color="grey30", 
                                        size=0.5, linetype="solid"),
        #panel.margin = unit(c(-0.5, -0.35, -0.5, -0.35), "cm"),
        #plot.margin = unit(c(0.1, 0.3, 0, 0.5), "cm")
  )
p4


## Combine the plots into a single layout

# Set up the page
sq.main <- nlevels(plotterBC[[1]])*4

require(patchwork)
confmatpl <- p1 +  p3 + p2 +p4 + 
  plot_layout(heights=c(sq.main, 4.7),widths=c(sq.main, 4.7))
confmatpl

### Importance plot ----
imppl <- data.table(Var=rownames(imps[[1]][[1]]))

for (i in 1:10){
  
  imppl <- cbind(imppl,as.data.table(imps[[i]][[1]])[,1])
  
}

setnames(imppl,c('Var','Imp1','Imp2','Imp3','Imp4','Imp5','Imp6','Imp7','Imp8','Imp9','Imp10'))

imppl[,
      c("Mean",'Sd','Se') := 
        .(rowMeans(.SD, na.rm = TRUE), 
          apply(.SD, 1, sd, na.rm = TRUE),
          apply(.SD, 1, plotrix::std.error, na.rm = TRUE)), 
      .SDcols = 2:11]

imppl[,Var:=factor(Var,levels=Var[order(Mean)])]

## Add plot to importance plot data to collator list
VIs <- rbind(VIs,data.table(Tax=tax,imppl[,.SD,.SDcols=c('Var','Mean','Se')]))

impplot <-  ggplot(imppl,(aes(x=Var,y=Mean))) +
  geom_bar(stat = 'identity',fill='#66afad', col='#172957',) +
  scale_x_discrete(labels=envlab[levels(imppl$Var)]) +
  geom_linerange(inherit.aes=FALSE,
                 aes(x=Var, ymin=Mean-Se, ymax=Mean+Se), 
                 colour='#172957', alpha=0.9, linewidth=1.3) +
  ylab(label = 'Mean decrease in Gini coefficient') +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5,hjust = 1, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_text(vjust=-4, size=12,colour="black"),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey90"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"),
        plot.margin = ggplot2::margin(0.5, 0.5, 1, 0.5, "cm"),)
impplot


## Response plot ----

# List for plot objects
cvRP.list <- list()


for (i in predselnf) {
  
  cvRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y,col=class)) +
    geom_smooth(method = 'loess',size=0.8,se=FALSE,span = 0.3,show.legend = FALSE) +
    facet_wrap(~ predvar,scales = "free_x", ncol=3,labeller = labeller(predvar=envlab)) +
    scale_colour_manual(values = classpal.r) +
    ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
    theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
          axis.title.x  = element_blank(),
          plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
          strip.background = element_rect(fill="grey90"),
          strip.text.x = element_text(size=12, face="bold"),
          panel.grid.major = element_line(colour="grey80"),
          panel.grid.minor = element_line(colour="grey80"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=12,colour="black"))
}

cvRP.list[['GearType']] <-  ggplot(plotdata2, aes(x=x,y=y,fill=class)) +
  geom_bar(stat = 'identity',position = "dodge",show.legend = FALSE) +
  scale_fill_manual(values = classpal.r) +
  facet_wrap(~ predvar,scales = "free_x", ncol=3,labeller = labeller(predvar=envlab)) +
  ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
  theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black",angle = 0),
        axis.title.x  = element_blank(),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey90"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"))

cvRP.list[['Legend']] <- as_ggplot(get_legend(ggplot(plotdata[plotdata$predvar== "SPM_MEAN",],aes(x=x,y=y,col=class)) +
                                                geom_smooth(method = 'loess') +
                                                scale_colour_manual(values = classpal.r,name='Response trait class',
                                                                    guide=guide_legend(title.position = 'top',)) +
                                                theme(legend.direction = 'horizontal',
                                                      legend.key = element_blank(),
                                                      axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
                                                      axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
                                                      axis.text.x  = element_text(vjust=0.5, size=12,colour="black",angle = 0),
                                                      axis.title.x  = element_blank(),
                                                      plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
                                                      strip.background = element_rect(fill="grey90"),
                                                      strip.text.x = element_text(size=12, face="bold"),
                                                      panel.grid.major = element_line(colour="grey80"),
                                                      panel.grid.minor = element_line(colour="grey80"),
                                                      panel.background = element_rect(fill="white"),
                                                      #legend.key = element_rect(fill = NA),
                                                      #legend.text = element_text(size=12,colour="black")
                                                )))


# Combine all plots to one canvas
cvRP <- wrap_plots(cvRP.list) +
  plot_layout(heights=1,widths=2)
cvRP


# Whole model performance table
asg.perf[, data.table(t(.SD), keep.rownames=TRUE),] %>%
  kbl('html',digits = 2,escape = FALSE, col.names = c('Statistic','Mean \u00B1 SD'),
      caption='Whole model performance') %>%
  kable_classic(full_width = F, position = "left",fixed_thead = T) %>%
  row_spec(0, bold = T)  %>%
  column_spec(1:2, width = "3cm") 

# Class-specific performance table
asg.cl.perf[-7, ] %>%
  kbl('html',digits = 2,escape = FALSE, align=c('l',rep('r', 4)),
      caption='Class-specific performance') %>%
  kable_classic(full_width = F, position = "left",fixed_thead = T) %>%
  row_spec(0, bold = T)  %>%
  column_spec(1:2, width = "1.5cm") %>%
  column_spec(3:5, width = "3cm")

### Make model predictions ----

# Drop unnecessary predictor layers 
dr <- names(s1)
dr <- dr[!dr %in% names(ffs[[1]]$forest$xlevels)]
s2 <- dropLayer(s1, dr)

# Set gear type as standard across prediction
GearType <- data.frame(GearType=factor('Scoop',levels = ffs[[1]]$forest$xlevels$GearType))

# Set up objects to save rasters into
cvpred <- NULL
cvpred.cps <- list()

# Set number of runs
nruns=10

# Predict to all model runs
for (i in 1:length(ffs)){
  
  rnn <-  paste0('Run',i)
  
  if (is.null(cvpred)){
    cvpred  <- stack(raster::predict(s2,ffs[[i]],const=GearType))
    names(cvpred) <- rnn
  } else {
    tmpl <- predict(s2,ffs[[i]],const=GearType)
    names(tmpl) <- rnn
    cvpred <- addLayer(cvpred,tmpl)
  }
  
  
  cvpred.cps[[rnn]] <- predict(s2,ffs[[i]],const=GearType,type='prob',index=1:numclass.r)
  
}



### Create a raster stack for spatial confidence results
ROutput <- stack()

### Calculate most frequent class and its frequency

# Most frequent class
MaxClass <- modal(cvpred,freq=FALSE)
ROutput <- addLayer(ROutput,MaxClass)
# Frequency of most frequent class (fraction of runs)
MaxClassF <- modal(cvpred,freq=TRUE)/nruns
ROutput <- addLayer(ROutput,MaxClassF)

### Calculate average probabilities for classes
classsums <- Reduce("+", cvpred.cps)
AvePclass <- classsums / nruns

### Find average probability of maximum frequency class
MaxClassAveProb <- stackSelect(AvePclass, MaxClass)
ROutput <- addLayer(ROutput,MaxClassAveProb)

### Calculate new layer for frequency x probability
CombConf <- MaxClassF * MaxClassAveProb
ROutput <- addLayer(ROutput,CombConf)

### Rename layers
names(ROutput) <- c("MaxClass","MaxClassF","MaxClassAveProb","CombConf")

### Plot layers 
plot(ROutput)

rm(MaxClass,MaxClassF,classsums,AvePclass,MaxClassAveProb,CombConf)

### Export Raster
writeRaster(ROutput, "ResponseTraitRasterBootOutput_Apr22.tif", format="GTiff",overwrite=T)
writeRaster(ROutput$MaxClass, "ResponseTraitMaxClass_Apr22.tif", format="GTiff",overwrite=T)
writeRaster(ROutput$MaxClassF, "ResponseTraitMaxClassFrequency_Apr22.tif", format="GTiff",overwrite=T)
writeRaster(ROutput$MaxClassAveProb, "ResponseTraitMaxClassAveProb_Apr22.tif", format="GTiff",overwrite=T)
writeRaster(ROutput$CombConf, "ResponseTraitConfidence_Apr22.tif", format="GTiff",overwrite=T)

### Save raster stack to R workspace
save(AvePclass,ROutput,file="ResponseTraitRasterBootOutput.RData")


# Summary statistics for the rasters
ROutput$MaxClass
ROutput$MaxClass <- as.factor(ROutput$MaxClass)
rat <- levels(ROutput$MaxClass)[[1]]
rat$Cluster <- levels(pshp$resptrait)
levels(ROutput$MaxClass) <- rat


# Plot rasters
extbb <- tmaptools::bb(pshp_sf,ext = 1.1,)
mpr <-  tm_shape(ROutput$MaxClass) + 
  tm_raster(title = 'Response trait cluster', #col='Cluster',
            style="fixed",
            palette=classpal.r,
            breaks = c(1:7),
            labels = levels(pshp$resptrait),
            legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = c("right","bottom"),
            legend.hist.height = 0.3,
            legend.hist.width = 0.7) +
  tm_shape(countries) +
  tm_grid(lines = FALSE,n.x = 4, n.y = 3) +
  tm_polygons() +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom"))
mpr

msd <-  tm_shape(ROutput$CombConf) + 
  tm_raster(title = 'Confidence',
            style="fixed",
            breaks = c(0,0.2,0.4,0.6,0.8,1),
            palette="RdYlGn",
            n = 5,
            legend.hist = FALSE) +
  tm_layout(legend.outside = TRUE) +
  tm_shape(countries) +
  tm_grid(lines = FALSE,n.x = 4, n.y = 3) +
  tm_polygons()
msd




##### 5. Effects traits model   ########################################################

### Repeat above with ####
tax='efftraits'

#_______________________________________________________________________________
#### PRODUCE SIDE BY SIDE PLOT FOR RESPONSE TRAITS ####

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/')

## Load libraries
library(raster)
library(leaflet)
## Load raster data layers
bathy <- raster("DATA/Rasters/Final/bathy3.tif")
pr = raster('Y:/C8381_OWEC_POSEIDON/Working_Area/SDM/OUTPUT/OneBenthicModelPredictionsApril22/ResponseTraitMaxClass_Apr22.tif')

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
prc = raster('Y:/C8381_OWEC_POSEIDON/Working_Area/SDM/OUTPUT/OneBenthicModelPredictionsApril22/ResponseTraitConfidence_Apr22.tif') #Anna's version

cols <- colorBin(palette = brewer.pal(n = 4, name = "Oranges"), bins = c(0,0.10,0.25,0.6), domain = c(0,0.10,0.25,0.6),na.color = "transparent")


plot(hill,col=grey(1:100/100),legend=FALSE)#,axes=F,box=FALSE
plot(pr, col=colours,add=TRUE,alpha=0.7,axes=F,box=FALSE,legend=FALSE)#no axes

#my.palette <- brewer.pal(n = 5, name = "Oranges")
my.palette <- c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603")

cuts=c(0,0.2,0.4,0.6,0.8,1.0) #set breaks,1.0


## SIDE BY SIDE PLOTS WITH a) and b) but on PLOT

png('C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTraitsMapping/OUTPUTS/FIGURE_2.jpg',width = 32,height = 15.80, units = "cm", res = 800,pointsize = 14)
# 2. Create the plot
line = 1
cex = 1.5
side = 3
adj=-0.15

par(mfrow=1:2)
#par(mar= c(2, 3.1, 2, 4)+ 0.2)
# margins are in the order: bottom, left, top, right
par(mar= c(4, 4, 2, 1)+ 0.2)
plot(pr, col=colours,add=F,alpha=1,axes=T,box=T,legend=F,xlab="Longitude", ylab="Latitude")
#legend(x = 8, y = 53.2, legend = c("2","1","3","4","5","6"), fill = c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59'),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE
legend(x = 7.9, y = 52.8, legend = c("2","1","3","4","5","6"), fill = c('#00E600','#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59'),cex = 1, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE

#
mtext("a)", side=side, line=line, cex=cex, adj=adj)
plot(prc,breaks=cuts,col=my.palette,add=F,alpha=1,axes=T,box=T,legend=F,xlab="Longitude", ylab="Latitude")
#legend(x = 4.8, y = 52.5, legend = c("0.0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1.0"), fill = c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),cex = 1.2, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE
legend(x = 5.1, y = 52.1, legend = c("0.0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1.0"), fill = c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),cex = 1, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE

mtext("b)", side=side, line=line, cex=cex, adj=adj)
dev.off()


#_______________________________________________________________________________
#### PRODUCE SIDE BY SIDE PLOT FOR EFFECTS TRAITS ####

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/')

## Load libraries
library(raster)

## Load raster data layers
bathy <- raster("DATA/Rasters/Final/bathy3.tif")
pr = raster('Y:/C8381_OWEC_POSEIDON/Working_Area/SDM/OUTPUT/OneBenthicModelPredictionsApril22/EffectsTraitMaxClass_Apr22.tif')

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
#prc = raster('Y:/C8210_OneBenthic/Working_Area/SDM/OUTPUT/EffectsTraitConfidence.tif') #Anna's version
prc = raster('Y:/C8381_OWEC_POSEIDON/Working_Area/SDM/OUTPUT/OneBenthicModelPredictionsApril22/EffectsTraitConfidence_Apr22.tif') #Anna's version
#cols <- colorBin(palette = brewer.pal(n = 4, name = "Oranges"), bins = c(0,0.10,0.25,0.6), domain = c(0,0.10,0.25,0.6),na.color = "transparent")

plot(hill,col=grey(1:100/100),legend=FALSE)#,axes=F,box=FALSE
plot(pr, col=colours,add=TRUE,alpha=0.7,axes=F,box=FALSE,legend=FALSE)#no axes

#my.palette <- brewer.pal(n = 4, name = "Oranges")
my.palette <- c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603")

cuts=c(0,0.2,0.4,0.6,0.8,1.0) #set breaks,1.0

## SIDE BY SIDE PLOTS WITH a) and b) but on PLOT

#png('OUTPUTS/EFFECTS_TRAITS/traits_effects_modelplot_confidenceplot.png',width = 30,height = 13.5,units = "cm", res = 600,pointsize = 12)
png('C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTraitsMapping/OUTPUTS/FIGURE_4.jpg',width = 32,height = 15.80,units = "cm", res = 800,pointsize = 14)

# 2. Create the plot
line = 1
cex = 1.5
side = 3
adj=-0.15

par(mfrow=1:2)
par(mar= c(4, 4, 2, 1)+ 0.2)
plot(pr, col=colours,add=F,alpha=1,axes=T,box=T,legend=F,xlab="Longitude", ylab="Latitude")#
legend(x = 7.9, y = 52.8, legend = c("6","3","1","2","4","5"), fill = c('#FF0000', '#c7eae5','#5ab4ac', '#d8b365', '#8c510a', '#f6e8c3'),cex = 1, inset = 0.9,bty = "n") # bty used to turn off legend border) EFFECTS
mtext("a)", side=side, line=line, cex=cex, adj=adj)
plot(prc,col=my.palette,add=F,alpha=1,axes=T,box=T,legend=F,xlab="Longitude", ylab="Latitude")
legend(x = 5.1, y = 52.1, legend = c("0.0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1.0"), fill = c("#FEEDDE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),cex = 1, inset = 0.9,bty = "n") # bty used to turn off legend border) RESPONSE
mtext("b)", side=side, line=line, cex=cex, adj=adj)

dev.off()
#_______________________________________________________________________________