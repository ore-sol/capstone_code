rm(list=ls())

## -------------- load in packages ## -------------- 

packages <- list('dplyr', 'terra','sf','vegan','ggplot2', 'tidyr','conflicted',
                 'landscapemetrics', 'BiocManager', 'phyloseq','reshape2')
lapply(packages,library, character.only=TRUE)

## -------------- take in csv files from output ## --------------

field_metadata <- read.csv('../Data/wheat_kits/sample_field_data1.csv')

#import and format colnames to match metadata
fungi2 <- read.csv('../Data/wheat_kits/2_UREN_WP1_ITS_COLLAPSED_99_FARMKITS_ONLY.csv')
colnames(fungi2) <- gsub("Sample_","", colnames(fungi2)) 
colnames(fungi2) <- gsub("hifi_reads_.*","hifi_reads",colnames(fungi2)) 
colnames(fungi2) <- gsub("..bc","--bc",colnames(fungi2), fixed = TRUE) 

new.fungi2 <- melt(fungi2)
#new.fungi2 <- subset(new.fungi2, select = -variable) %>% unique() 
new.fungi2 <- new.fungi2[!is.na(new.fungi2$family),] 
new.fungi2 <- new.fungi2[new.fungi2$family !="",] 


#import and format colnames to match metadata
fungi <- read.csv('../Data/wheat_kits/UREN_WP1_ITS_COLLAPSED_99_FARMKITS_ONLY.csv')
colnames(fungi) <- gsub("Sample_","", colnames(fungi)) 
colnames(fungi) <- gsub("hifi_reads_.*","hifi_reads",colnames(fungi)) 
colnames(fungi) <- gsub("..bc","--bc",colnames(fungi), fixed = TRUE) 

new.fungi <- melt(fungi)
#new.fungi <- subset(new.fungi, select = -variable) %>% unique() 
new.fungi <- new.fungi[!is.na(new.fungi$family),] 
new.fungi <- new.fungi[new.fungi$family !="",] 

#import and format colnames to match metadata
bacteria <- read.csv('../Data/wheat_kits/UREN_16S_Collapsed_99_updated.csv')
colnames(bacteria) <- gsub("Sample_","", colnames(bacteria)) 
colnames(bacteria) <- gsub("hifi_reads_.*","hifi_reads",colnames(bacteria)) 
colnames(bacteria) <- gsub("..bc","--bc",colnames(bacteria), fixed = TRUE) 

new.bacteria <- bacteria[!grep("*hifi_reads*",colnames(bacteria))] %>% melt()
new.bacteria <- melt(bacteria)
new.bacteria <- subset(new.bacteria, select = -variable) %>% unique() 
new.bacteria <- new.bacteria[!is.na(new.bacteria$family),] 
new.bacteria <- new.bacteria[new.bacteria$family !="",] 

#read what columns are there and relevant for: , generally looking through the data 
colnames(fungi)
colnames(fungi2)
colnames(bacteria)
colnames(field_metadata)

#fungi[fungi$ESV=='ESV_7474',] %>% View()

## -------------- separate files into by batch 1 or batch 2 ...? ## -------------- 

## -------------- make alternate dfs where you have 1) filtered the fungi/bacteria df by origin status and then 2) grab those identications
# -> use metadata file for grouping decision ... 

####take all grouping: 
takeall <- field_metadata[!field_metadata$take_all_seen =="unknown" & !is.na(field_metadata$take_all_seen),]

#wrangle data 
takeall<- takeall[c('SampleID','take_all_seen')] #now you can join to fungi / bacterial df


## ultimately i think an if statement might just be faster & more reprod. than all this ...




## -------------- rank soils by microorg composition (relative) ## -------------- 

## -------------- feed ranking data into random forest machine learning model ## -------------- 

## -------------- make map of where plots/samples were collected from ## -------------- 
sites <- field_metadata[!is.na(field_metadata$gps_latitude) & !is.na(field_metadata$gps_longitude),] %>%
  
  st_as_sf(coords=c('gps_longitude','gps_latitude'), crs=4326)

nrow(sites)
nrow(field_metadata)
# a map.tif of the area.... 
landcover <- terra::rast("../Data/uk_srtm_shade.tif") 
print(landcover) # map details, interesting to note: each numerical value is a land cover type

# Get the bounding box and convert to a polygon
sites_region <- st_as_sfc(st_bbox(sites)) #makes a box capturing all of your points
# Buffer the region by 0.1 degrees to get coverage around the landscape.
# This is roughly 10 kilometres which is a wide enough border to allow
# us to test the scale of effects.
sites_region <- st_buffer(sites_region, 0.01) #adds some buffer area to box made on line

sites_utm23S <- st_transform(sites, 32723) # converted vector to another coord. system


plot(landcover_utm23S) # <-- !! FIND OUT WHY THIS DOESNT WORK
plot(st_geometry(sites_utm23S), add=TRUE)
