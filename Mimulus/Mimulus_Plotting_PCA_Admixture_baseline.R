#############################
# R plotting
############################
# - List of files downloaded from Cedar
# - Package installations
# - libraries
# - Load in all the data 

# PCA analysis

# Join all the data into population file and sample file

# Admixture
# - Determine Admixture groups order
# - Set groups colours
# - Plot PCA images with colours
# - Plot map of sites with dominant colours
# - Make Admixture barplot (total, Ancient, Ellemere)
# - Make Admixture maps (samples and population scales)
# - Kluane correlation and barplot

# PopStats
# - Icetime, Diversity and Temperature correlations
# - Plots of Pi and Tajima's D per population
# - Maps of all population avgeraged stats

# Table of PopStats for paper

# Site specific Admixture plots

#####################################
# download files from Cedar to computer

# Admixture
# .4.Q

# # PopStats
# Fst_mean.txt
# Fst_weighted.txt
# Fst_count.txt
# Fst_sites.txt
# population_TajimaD.txt
# summed_site_Pi.txt
# summed_wind_Pi.txt
# Het_data_by_pop.txt

# # PCA GDS file
# Cassiope_noMER_r10i.recode.gds


######################################

# # install packages
# install.packages("devtools") 
# devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
# devtools::install_github("celphin/Equitable")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# 
# install.packages(c("tidyverse","dplyr","vctrs", "tidyr","reshape","geosphere","ade4","ape","scatterpie","data.table","maps","readr","stringr","ggplot2"))

#----------------------------------
# libraries

library(Equitable)
library(PlotSvalbard)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape)
library(geosphere)
library(ade4)
library(ape)
library(scatterpie)
library(data.table)
library(maps)
library(readr) 
library(stringr)
library(SNPRelate)
library(ggplot2)

##############################
# Load in the data

setwd("~/GitHub/Population_genomics_Cassiope/Mimulus/Baseline")

#Load the gds file - for PCA
genofile <- snpgdsOpen("./Figures_data/Mimulus_filtered_baseline.recode.gds")

# load list of samples in filtered vcf
samples_list <- read.table("./Figures_data/Mimulus_filtered_baseline.imiss", header = TRUE)

# load in Admixture data - code setup for 5 populations
Admix_tbl=read.table("./Figures_data/Admixture/Mimulus_filtered_baseline.4.Q")

# load detailed information about all 371 individuals
sample_information  <- read.csv("./Figures_data/M_caridnalis_genomics_meta_data.csv", header = TRUE)

# Population lat long
pop_lat_long <- read.csv("./Figures_data/56_pop_lat.csv", header = TRUE)



# PopStats
#Het_data <- read.table("./Figures_data/PopStats/Het_data_by_pop.txt", header = TRUE)
#Pi_wind_data <- read.table("./Figures_data/PopStats/summed_wind_Pi.txt", header=TRUE, sep=" ")
#Pi_site_data <- read.table("./Figures_data/PopStats/summed_site_Pi.txt", header=TRUE, sep=" ")
#Tajima_data <- read.table("./Figures_data/PopStats/population_TajimaD.txt", header=TRUE, sep=" ")


#################################
# PCA
# https://owensgl.github.io/biol525D/Topic_8-9/pca.html

# create GDS file on server
# snpgdsVCF2GDS("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.vcf.gz",
# "/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.gds",
# method="biallelic.only")

#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

#Make a list of sites we're keeping.
snpset.id <- unlist(snpset_pruned)

#Run the PCA
pca0 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, autosome.only = F)
pca5 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, maf=0.05,  autosome.only = F)

###################################
pca <- pca5

#Here's the percent variance explained for each eigenvector
pc.percent <- pca$varprop*100
round(pc.percent, 2)
# 19.99 18.55 15.99 13.10 11.25  6.45  4.42  2.81  2.58  1.11  1.05  0.84  0.67  0.42  0.39  0.26 

#Make a dataframe of your PCA results
PCA_tab <- data.frame(sample = pca$sample.id,
                      PC1 = pca$eigenvect[,1],    # the first eigenvector
                      PC2 = pca$eigenvect[,2],    # the second eigenvector
                      PC3 = pca$eigenvect[,3],    # the first eigenvector
                      PC4 = pca$eigenvect[,4],    # the second eigenvector
                      PC5 = pca$eigenvect[,5],    # the first eigenvector
                      PC6 = pca$eigenvect[,6],    # the second eigenvector
                      stringsAsFactors = FALSE)

str_split(as.character(PCA_tab$sample), "_")

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(PCA_tab$sample), "_"))))))
PCA_tab_data <- cbind(PCA_tab, q[,1])

colnames(PCA_tab_data) <- c("ID_code", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6","Pop")

#####################################
# Join all the data

# Sample specific information

# samples_list 
# Admix_tbl
# sample_information 
# Het_data
# PCA_tab_data

# Join Samples information
samples_list$ID_code <- as.factor(samples_list$INDV)
ID_code = samples_list$ID_code

Admix_tbl1 <- cbind(ID_code, Admix_tbl)
All_samples_data <- left_join(Admix_tbl1, sample_information, by="ID_code")
#All_samples_data <- left_join(All_samples_data, Het_data, by="ID_code")
All_samples_data <- left_join(All_samples_data, samples_list, by="ID_code")
All_samples_data <- left_join(All_samples_data, PCA_tab_data, by="ID_code")

#------------------------------
# join all Population information

# population_information  
# Pi_wind_data 
# Pi_site_data 
# Tajima_data 
# 
# colnames(Pi_site_data) <- c("Pop" ,   "SumSitePi" , "Site_SD_Pi",  "Site_var_Pi")
# colnames(Pi_wind_data) <- c("Pop" ,   "SumWindPi",  "Wind_SD_Pi",  "Wind_var_Pi")
# colnames(Tajima_data ) <- c("Pop" ,   "TajimaAvg",  "TajimaSD",  "TajimaVar")
# population_information$Pop[which(population_information$Pop=="25")] <- "025"
# 
# All_pop_data <- left_join(population_information, Tajima_data, by="Pop")
# All_pop_data <- left_join(All_pop_data, Pi_site_data, by="Pop")
# All_pop_data <- left_join(All_pop_data, Pi_wind_data, by="Pop")

#------------------------------
# average the sample information by population
Pop_avg <- All_samples_data  %>%
  group_by(Site) %>%
  dplyr::summarise(V1 = mean(V1, na.rm=TRUE), 
                   V2 = mean(V2, na.rm=TRUE), 
                   V3 = mean(V3, na.rm=TRUE), 
                   V4 = mean(V4, na.rm=TRUE), 
                   #V5 = mean(V5, na.rm=TRUE), 
                   # V6 = mean(V6, na.rm=TRUE), 
                   # V7 = mean(V7, na.rm=TRUE), 
                   # V8 = mean(V8, na.rm=TRUE), 
                   # V9 = mean(V9, na.rm=TRUE), 
                   # V10 = mean(V10, na.rm=TRUE),
                   # Lat = mean(Lat, na.rm=TRUE),
                   # Long = mean(Long, na.rm=TRUE),
                   # Lat_shift = mean(Lat_shift, na.rm=TRUE),
                   # Long_shift = mean(Long_shift, na.rm=TRUE),
                   # Avg.Temp = mean(Avg.Temp, na.rm=TRUE),
                   # Leaf_weight = mean(Leaf_weight, na.rm=TRUE),
                   # Ice_time = mean(Ice_time, na.rm=TRUE),
                   # O.HOM. = mean(O.HOM., na.rm=TRUE),
                   # E.HOM. = mean(E.HOM., na.rm=TRUE),
                   # N_SITES = mean(N_SITES, na.rm=TRUE),
                   # FIS = mean(F, na.rm=TRUE),
                   N_MISS = mean(N_MISS, na.rm=TRUE),
                   F_MISS = mean(F_MISS, na.rm=TRUE),
                   PC1 = mean(PC1, na.rm=TRUE),
                   PC2 = mean(PC2, na.rm=TRUE),
                   PC3 = mean(PC3, na.rm=TRUE),
                   PC4 = mean(PC4, na.rm=TRUE),
                   PC5 = mean(PC5, na.rm=TRUE),
                   PC6 = mean(PC6, na.rm=TRUE), 
                   Group = list(Geo_Region),
                   Paper_ID = Paper_ID
  )

Pop_avg$Pop <- Pop_avg$Site
All_pop_data <- left_join(pop_lat_long, Pop_avg, by="Paper_ID")

All_pop_data <- All_pop_data[-which(duplicated(All_pop_data)),] 

############################################
# calculate the max Admix group for each location
Admix_groups <- gather(All_pop_data, Group, Amount, V1:V4, factor_key=TRUE)

maxAdmixgroup <- Admix_groups  %>%
  group_by(Pop) %>%
  dplyr::summarise(Group = Group[which(Amount==max(Amount, na.rm=TRUE))])

Admix_groups$Mix <- paste0(Admix_groups$Group, Admix_groups$Pop)
maxAdmixgroup$Mix <- paste0(maxAdmixgroup$Group, maxAdmixgroup$Pop)

Final_AdmixGroups <- left_join(maxAdmixgroup, Admix_groups, by="Mix")


Max_Admix_Group <- c("V1", "V2", "V3", "V4")
Admix_Location <- c("North", "Center1", "Center2", "South")
Group_location <- as.data.frame(cbind(Max_Admix_Group, Admix_Location))

All_pop_data <- left_join(All_pop_data, Final_AdmixGroups, by="Site")

All_pop_data$Max_Admix_Group <- All_pop_data$Group.y
All_pop_data <- left_join(All_pop_data, Group_location, by="Max_Admix_Group")

#4 groups
Groups_summary_admix <- All_pop_data %>%
  group_by(Admix_Location) %>%
  dplyr::summarise(Locations = list(Pop.x), 
                   #Lat = mean(Lat.x, na.rm=TRUE),
                   #Long = mean(Long.x, na.rm=TRUE),
                   V1 = mean(V1, na.rm=TRUE), 
                   V2 = mean(V2, na.rm=TRUE), 
                   V3 = mean(V3, na.rm=TRUE), 
                   V4 = mean(V4, na.rm=TRUE))


###################################
# Set up colour scheme 

map_colours_5g <- c("deepskyblue", "yellow",  "green", "red")

Max_Admix_Group <- c("V1", "V2","V3", "V4")

map_col <- cbind(Max_Admix_Group, map_colours_5g)
map_col <- as.data.frame(map_col)

All_pop_data <- left_join(All_pop_data, map_col, by="Max_Admix_Group")

# join map_colours with individual sample data
Map_col_pop <- as.data.frame(cbind(All_pop_data$Pop, All_pop_data$map_colours_5g, All_pop_data$Max_Admix_Group))

colnames(Map_col_pop) <- c("Pop", "map_colours_5g", "Max_Admix_Group")

All_samples_data$Pop <- All_samples_data$Site

All_samples_data <- left_join(All_samples_data, Map_col_pop, by="Pop")

All_samples_data$map_colours_5g <- as.factor(All_samples_data$map_colours_5g)

#####################################

#Plot a PCA image

# PCA = 27.00  5.62  2.63  1.63  1.50  1.26  1.02  0.93  0.93  0.85  0.82  0.77  0.73  0.73  0.68  0.67
All_samples_data$Site <- as.factor(All_samples_data$Site)
jpeg("./Figures_data/Plots/PCA_PC1_PC2_maf1.jpg", width = 3000, height = 2700)
All_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = Max_Admix_Group), size=15)  +
  geom_text(aes(label = ID_code), color = "black", fontface = 2, size = 25.4/72.27*15)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC2", x = "PC1")
dev.off()

jpeg("./Figures_data/Plots/PCA_PC2_PC3_maf1.jpg", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC2,y=PC3)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC3 (2.63%)", x = "PC2 (5.62%)")
dev.off()

jpeg("./Figures_data/Plots/PCA_PC3_PC4_maf1.jpg", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC3,y=PC4)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC4 (1.63%)", x = "PC3 (2.63%)")
dev.off()

jpeg("./Figures_data/Plots/PCA_PC5_PC6_maf1.jpg", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC5,y=PC6)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC6 (1.26%)", x = "PC5 (1.50%)")
dev.off()

############################################

#Make Admixture barplot

#Excel: =RIGHT(A2,LEN(A2) - SEARCH("_", A2, SEARCH("_", A2) + 1))

All_lat_long_data <-  All_pop_data[,c(1:5)]
All_samples_data$Site <- All_samples_data$Pop

All_samples_data_map <- left_join(All_samples_data, All_lat_long_data, by="Site")

mergedAdmixTable <- All_samples_data_map[,c("Site", "ID_code","Latitude.x", "Longitude.x","V1", "V2","V3", "V4")]

rownames(mergedAdmixTable) <- mergedAdmixTable$ID_code

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

#----------------------------------
#4 Groups

ordered4 = mergedAdmixTable[order(mergedAdmixTable$Geo_Region),]

ordered0 = mergedAdmixTable[order(mergedAdmixTable$Site),]
ordered1 = ordered0[order(ordered0$V3),]
ordered2 = ordered1[order(ordered1$V2),]
ordered3 = ordered2[order(ordered2$V1),]
ordered4 = ordered3[order(ordered3$V4),]

ordered4 = mergedAdmixTable[order(mergedAdmixTable$Latitude.x),]


ordered4$Site <-  as.factor(ordered4$Site)

#map_colours_5g <- c("deepskyblue", "yellow",  "green", "red")
map_colours_5g <- c("green",  "yellow", "deepskyblue", "red")

jpeg("./Figures_data/Plots/Site_Admix_4_bar.jpg", width = 6000, height = 2000)
# bottom, left, top, and right
par(mar=c(30,10,4,4))
barplot(t(as.matrix(ordered4[,c(5:8)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(ordered4$ID_code), las=2, cex.names=1, cex.axis=6)
dev.off()


#########################################
# Map of samples
# https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html

# https://ggplot2-book.org/maps.html
# https://jakob.schwalb-willmann.de/basemaps/
All_pop_data  <- transform_coord(All_pop_data, lon = "Longitude.x", lat = "Latitude.x", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3857")
All_pop_data$Pop <- as.factor(All_pop_data$Pop)
#All_pop_data <- All_pop_data[,-c(21:22)]

library(basemaps)
library(mapedit)
library(ggmap)
library(ggnewscale)
# view all available maps
#get_maptypes()

# use draw_ext() to interactively draw an extent yourself
# ext <- draw_ext()

# load and return basemap map as class of choice, e.g. as image using magick:
#basemap_magick(ext, map_service = "osm", map_type = "topographic")

# set defaults for the basemap
#set_defaults(map_service = "osm_stamen", map_type = "toner")
#set_defaults(map_service = "osm_stamen", map_type = "terrain_bg")
#set_defaults(map_service = "osm", map_type = "topographic")
set_defaults(map_service = "osm_stamen", map_type = "terrain")

#All_pop_data <- left_join(All_pop_data, Map_col_pop, by="Pop")

jpeg("./Figures_data/Plots/Mimulus_samples_map_colours.jpg",width = 2700, height = 3300)
ggplot() + 
  basemap_gglayer(ext) +
  scale_fill_identity() + 
  coord_sf() +
  geom_point(data = All_pop_data, aes(x = lon.utm,  y = lat.utm, col = "Population"), color="yellow", size = 20) +
  geom_text(data = All_pop_data, aes(x = lon.utm, y = lat.utm, label = Site), color = "black", fontface = 2, size = 25.4/72.27*30)
dev.off()

#https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/

jpeg("./Figures_data/Plots/Admix_map4_take1_lat_long.jpg", width = 900, height = 1500)
ggplot(data = All_pop_data, aes(x = lon.utm, y = lat.utm)) +
  #basemap_gglayer(ext) +
  #coord_sf() +
  #scale_fill_identity() +
  #new_scale_colour() +
  geom_scatterpie(aes(x = lon.utm, y = lat.utm, group=Site, r=20000), data = All_pop_data, cols=colnames(All_pop_data[,c(6:9)]), size = 0.1)+
  scale_fill_manual(values=map_colours_5g) +
  geom_text(data = All_pop_data, aes(x = lon.utm, y = lat.utm, label = Site), color = "black", fontface = 2, size = 25.4/72.27*10) 
dev.off()

#cannot get the two scales for the map and the pie charts to work

*************************
  
##############################
# Write out samples file

write.table(All_samples_data, file = "./Figures_data/All_samples_data.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#####################################
# Admixture map - shifted to see all groups

#library(devtools)
#devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
#library(PlotSvalbard)
# https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html
#library(scatterpie)

# shifted coordinates
All_samples_data <- transform_coord(All_samples_data, lon = "Long_shift", lat = "Lat_shift", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")
All_samples_data$lon_shift.utm <- All_samples_data$lon.utm 
All_samples_data$lat_shift.utm <- All_samples_data$lat.utm 
All_samples_data <- subset(All_samples_data, select = -c(lon.utm, lat.utm))

# regular coordinates
All_samples_data <- transform_coord(All_samples_data, lon = "Long", lat = "Lat", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

# make dataset of just what needs to be plotted
shift_Admix_lat_long <- subset(All_samples_data, select = c(ID_code, lon.utm, lat.utm, lon_shift.utm, lat_shift.utm, Pop, V1, V2, V3, V4, V5))

Pop <- c("V1", "V2", "V3", "V4", "V5")

jpeg("./Figures_data/Plots/Admix_map5_take1_shiftedlat_long.jpg", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon_shift.utm, y = lat_shift.utm, group = ID_code, r=200000), data = shift_Admix_lat_long, cols = Pop, size = 0.9) +
  geom_point(aes(x = lon_shift.utm, y = lat_shift.utm), data = shift_Admix_lat_long, col="white", size=17) +
  geom_point(aes(x = lon.utm, y = lat.utm), data = shift_Admix_lat_long, col="green", size=6) +
  scale_fill_manual(values=map_colours_5g)+
  geom_text(data = shift_Admix_lat_long, aes(x = lon_shift.utm, y = lat_shift.utm, label = Pop), color = "black", fontface = 2, size = 25.4/72.27*20)
dev.off()

#---------------------------------
# Averaged admixture map - not plotting each individual

Pop <- c("V1", "V2", "V3", "V4", "V5")

# map showing population structure makeup
jpeg("./Figures_data/Plots/Admix_map5g_take1_avg.jpg", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon.utm, y = lat.utm, group = Pop, r=200000), data = All_pop_data, cols = Pop, size = 0.5) +
  scale_fill_manual(values=map_colours_5g)
dev.off()

#----------------------------
All_pop_data <- subset(All_pop_data, select = -c(lon.utm, lat.utm))

All_pop_data <- transform_coord(All_pop_data, lon = "Long_shift", lat = "Lat_shift", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")
All_pop_data$lon_shift.utm <- All_pop_data$lon.utm 
All_pop_data$lat_shift.utm <- All_pop_data$lat.utm 
All_pop_data <- subset(All_pop_data, select = -c(lon.utm, lat.utm))

All_pop_data <- transform_coord(All_pop_data, lon = "Long.x", lat = "Lat.x", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

jpeg("./Figures_data/Plots/Admix_map5_take1_shiftedlat_long_avg.jpg", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon_shift.utm, y = lat_shift.utm, group = Pop, r=200000), data = All_pop_data, cols = Pop, size = 0.9) +
  geom_point(aes(x = lon_shift.utm, y = lat_shift.utm), data = All_pop_data, col="white", size=17) +
  geom_point(aes(x = lon.utm, y = lat.utm), data = All_pop_data, col="green", size=6) +
  scale_fill_manual(values=map_colours_5g)+
  geom_text(data = All_pop_data[-which(All_pop_data$Pop=="AlexOld"), ], aes(x = lon_shift.utm, y = lat_shift.utm, label = Pop), color = "black", fontface = 2, size = 25.4/72.27*20)
dev.off()

######################################

# Kluane
#check correlation with %BC and elevation

All_samples_data
Kluane <- All_samples_data[which(All_samples_data$Pop.x=="KL"|All_samples_data$Pop.x=="PC"),]

# Saximontana is V5 for Take 1

mod1 = lm(V3~Elevation, data = Kluane)
modsum = summary(mod1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.292e+00  1.440e-01   8.971 7.42e-08 ***
#   Elevation   -6.859e-04  9.243e-05  -7.421 1.00e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07808 on 17 degrees of freedom
# Multiple R-squared:  0.7641,	Adjusted R-squared:  0.7502 
# F-statistic: 55.07 on 1 and 17 DF,  p-value: 9.995e-07

jpeg("./Figures_data/Plots/Elevation_BC_5.jpg", width = 3000, height = 1400)
par(mar=c(20,20,4,4))
plot(Kluane$Elevation, Kluane$V3, pch=20, cex = 10, mgp=c(10,5,0), cex.lab=5, cex.axis=5, xlab="Elevation(m)", ylab="BC Admixture Proportion")
abline(mod1)
dev.off()


orderedKluane = Kluane[order(Kluane$Elevation),]

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

jpeg("./Figures_data/Plots/Admix_Kluane_bar5.jpg", width = 1000, height = 707)
barplot(t(as.matrix(orderedKluane[,c(3:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(orderedKluane$ID_code), las=2, cex.names=1.7)
dev.off()

#################################################
# Climate, ice and nucleotide diversity
# Icetime.x
# tavglong
# SumWindPi

# check site and window pi correlate
jpeg("./Figures_data/Plots/SitevsWIndPi.jpg", width = 1000, height = 707)
plot(All_pop_data$SumWindPi ~ All_pop_data$SumSitePi, pch=20, xlab="Sites Pi", ylab="Window Pi")
dev.off()

# Models
mod3 = lm(SumWindPi ~ tavglong, data = All_pop_data)
summary(mod3)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -14.2604  -1.5667   0.1728   1.9179  11.5762 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  16.8080     1.7931   9.374 5.95e-11 ***
#   tavglong      0.5029     0.2390   2.104   0.0428 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.205 on 34 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.1152,	Adjusted R-squared:  0.08918 
# F-statistic: 4.427 on 1 and 34 DF,  p-value: 0.04285


mod4 = lm(Ice_time.x ~ tavglong, data = All_pop_data)
summary(mod4)

# Residuals:
#   Min     1Q Median     3Q    Max 
# -9.730 -2.978 -1.069  1.762  8.380 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  10.9930     1.5866   6.928 3.07e-08 ***
#   tavglong      0.1624     0.1863   0.872    0.389    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.565 on 38 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.01962,	Adjusted R-squared:  -0.006184 
# F-statistic: 0.7603 on 1 and 38 DF,  p-value: 0.3887


mod5 = lm(SumWindPi ~ Ice_time.x, data = All_pop_data)
summary(mod5)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11.3903  -1.5145   0.0607   1.0588  14.5881 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  15.7977     1.8528   8.526 3.66e-10 ***
#   Ice_time.x    0.3520     0.1417   2.484   0.0178 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.082 on 36 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.1463,	Adjusted R-squared:  0.1226 
# F-statistic: 6.171 on 1 and 36 DF,  p-value: 0.01777


# Plots
jpeg("./Figures_data/Plots/WindPivsTemp.jpg", width = 1000, height = 707)
plot(All_pop_data$SumWindPi ~ All_pop_data$tavglong, col="lightblue", pch=19, cex=2, xlab="Average Temperature (degrees Celsius)", ylab="Sites Pi")
abline(mod3, col="red", lwd=3)
text(SumWindPi ~ tavglong, labels=Pop, data=All_pop_data, cex=0.9, font=2)
dev.off()

jpeg("./Figures_data/Plots/Ice_temp.jpg", width = 1000, height = 707)
plot(All_pop_data$Ice_time.x ~ All_pop_data$tavglong, col="lightblue", pch=19, cex=2, xlab="Average Temperature (degrees Celsius)", ylab="Icetime")
abline(mod4, col="red", lwd=3)
text(Ice_time.x ~ tavglong, labels=Pop,data=All_pop_data, cex=0.9, font=2)
dev.off()

jpeg("./Figures_data/Plots/IcetimesWIndPi.jpg", width = 1000, height = 707)
plot(All_pop_data$SumWindPi ~ All_pop_data$Ice_time.x, col="lightblue", pch=19, cex=2, xlab="Ice retreat time (thousands of years)", ylab="Sites Pi")
abline(mod5, col="red", lwd=3)
text(SumWindPi ~ Ice_time.x, labels=Pop, data=All_pop_data, cex=0.9, font=2)
dev.off()

##################
# make plots of Pi, Tajimas D

# order population averages by longitude
All_pop_data_ord <- All_pop_data[order(All_pop_data$Long.x),]

# remove Sverdrup New
All_pop_data_ord_sub <- All_pop_data_ord[-which(All_pop_data_ord$Pop=="SVN"),]

# remove NA values
All_pop_data_ord_sub <- All_pop_data_ord_sub[-which(is.na(All_pop_data_ord_sub$SumWindPi)),]

# make Pop a factor
All_pop_data_ord_sub$Pop <- as.factor(All_pop_data_ord_sub$Pop)

# plot Nucleotide diversity
jpeg("./Figures_data/Plots/AvgPi_noSverdrup.jpg", width = 2200, height = 707)
plot(All_pop_data_ord_sub$Pop, All_pop_data_ord_sub$SumWindPi)
dev.off()

jpeg("./Figures_data/Plots/AvgTajima_noSverdrup.jpg", width = 1700, height = 707)
plot(All_pop_data_ord_sub$Pop, All_pop_data_ord_sub$TajimaAvg, cex.axis=2, cex.lab=1.5, las=2)
dev.off()

##################
# make barplots of heterozygotsity with per individual data

# All_samples_data$F, All_samples_data$O.HOM. , All_samples_data$E.HOM.

All_samples_data$Pop <- as.factor(All_samples_data$Pop)

# jpeg("./Figures_data/Plots/FIS.jpg", width = 2200, height = 707)
# plot(All_samples_data$Pop, All_samples_data$F, col="blue")
# dev.off()

jpeg("./Figures_data/Plots/Obser_Homozygosity.jpg", width = 2200, height = 707)
plot(All_samples_data$Pop, All_samples_data$O.HOM., col="blue")
dev.off()

jpeg("./Figures_data/Plots/Exp_Homozygosity.jpg", width = 2200, height = 707)
plot(All_samples_data$Pop, All_samples_data$E.HOM., col="blue")
dev.off()


#-------
# order by longitude
Het_data <-  as.data.frame(cbind(All_samples_data$Pop.x, All_samples_data$F, All_samples_data$Long, All_samples_data$map_colours_5g))
colnames(Het_data) <- c("Pop", "F", "Long", "map_colours_5g")
Het_data$F <- as.numeric(Het_data$F)
Het_data$Long <- as.numeric(Het_data$Long)

Het_data_ord <- Het_data[order(Het_data$Long),]

Het_data_ord$Pop <- as.factor(Het_data_ord$Pop)
levels(Het_data_ord$Pop)
Het_data_ord$Pop <- factor(Het_data_ord$Pop, levels = c("GEN", "ATQ", "BARD", "MAT", "MNT", "IMN", "SAG" , "DEN", "MIL",  "QHI",  "KL", "PC",  "PEA", 
                                                        "KUQ",   "YAM", "AXE", "EUR",  "FOS", "GF", "SVO" , "BY" , "AlexNew", "AlexOld" ,  "HAZ",  "Iq",
                                                        "CR", "Kik", "DLG", "DQG",  "025", "IG", "ZAC", "LON", "PET", "LAJ", "SW", "SAM", "YED"))

Het_cols <-  c("purple", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue" , "deepskyblue", "deepskyblue",  "deepskyblue",  "deepskyblue", "deepskyblue",  "deepskyblue", 
               "deepskyblue",   "deepskyblue", "deepskyblue", "yellow",  "yellow", "yellow", "yellow" , "yellow" , "yellow", "yellow" ,  "yellow",  "orange",
               "orange", "orange", "orange", "orange",  "orange", "orange", "yellow", "yellow", "yellow", "yellow", "yellow", "red", "red")

# jpeg("./Figures_data/Plots/FIS_ordered.jpg", width = 2200, height = 707)
# plot(Het_data_ord$Pop, Het_data_ord$F, col=Het_cols)
# dev.off()


# par(mar = c(bottom, left, top, right)) 
jpeg("./Figures_data/Plots/Fis_col_box.jpg", width = 4000, height = 1700)
par(mar= c(20,15,4,4))
plot(F ~ Pop,  data = Het_data_ord, xlab="", ylab="", cex=3, cex.lab=6, cex.axis=5, las=2, col = Het_cols)
stripchart(F ~ Pop, vertical = TRUE, method = "jitter", pch = 16,las = 2, col = "blue", cex = 1.5, data = Het_data_ord, add=TRUE)
lines(Het_data_ord$Pop, y=rep(0, length(Het_data_ord$Pop)), lwd=3)
dev.off()



