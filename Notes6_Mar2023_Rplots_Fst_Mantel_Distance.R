#####################################
# Entire Cassiope analysis - part 6
# March 2023
# FST and Mantel test
###########################
# File contains: 

 # - FST download
 # - Distance between populaitons calculation
 # - Plotting FST for populations and Admix results
 # - Mantel Test
 
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

setwd("~/GitHub/Rieseberg_Lab/Population_genomics_Cassiope")

# Population specific data
All_pop_data <- read.table("./Figures_data/All_pop_data.txt",  sep = "\t", dec = ".", header = TRUE)

# load Admixture groups FST values
Fst_Groups <- read.table("./Figures_data/Admixture/Take1/Fst5.txt", header = TRUE)

# FST by sampling location
count <- read.table("./Figures_data/PopStats/Fst_count.txt", header = FALSE)
sites <- read.table("./Figures_data/PopStats/Fst_sites.txt", header = FALSE)
mean <- read.table("./Figures_data/PopStats/Fst_mean.txt", header = FALSE)
weighted <- read.table("./Figures_data/PopStats/Fst_weighted.txt", header = FALSE)

length(unique(sites$V2))

mean_weight<- cbind(sites, count, mean, weighted)
Fst_data <- mean_weight[, c(2, 6, 17, 24)]


#################################
# Determine distance between populations

Admix_group_dist <- as.data.frame(cbind(All_pop_data$Pop, All_pop_data$Lat.x, All_pop_data$Long.x))

dt <- expand.grid.df(Admix_group_dist,Admix_group_dist)
colnames(dt) <- c("Pop", "Lat", "Long", "Pop_dest", "Lat_dest","Long_dest")

dt$Lat_dest <- as.numeric(dt$Lat_dest)
dt$Long_dest <- as.numeric(dt$Long_dest)
dt$Lat <- as.numeric(dt$Lat)
dt$Long <- as.numeric(dt$Long)
dt$dist <- (distGeo(matrix(c(dt$Long, dt$Lat), ncol = 2), matrix(c(dt$Long_dest, dt$Lat_dest), ncol = 2)))/1000

dt_sub <-as.data.frame(cbind(as.character(dt$Pop), as.character(dt$Pop_dest), dt$dist))

group_distances <- spread(dt_sub, V2, V3)

write.table(group_distances, file = paste0("./Figures_data/group_distances.txt"), quote = FALSE, row.names=FALSE, col.names=TRUE)

#-----------------------
# get the longitudinal distance 

dt$Lat <- as.numeric(as.character(dt$Lat))
dt$Long <- as.numeric(as.character(dt$Long))
dt$Lat_dest <- as.numeric(as.character(dt$Lat_dest))
dt$Long_dest <- as.numeric(as.character(dt$Long_dest))

dt$dist2 <- sqrt(((dt$Lat - dt$Lat_dest)^2)+((dt$Long - dt$Long_dest)^2))
dt_sub2 <-as.data.frame(cbind(as.character(dt$Pop), as.character(dt$Pop_dest), dt$dist2))
group_distances2 <- spread(dt_sub2, V2, V3)


########################
# FST plotting

colnames(Fst_data) <- c("Sites","Count","MeanFst", "WeightedFst")

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(Fst_data$Sites), "_"))))))
colnames(q) <- c("FST", "Pop1", "Pop2")
Fst_data_pop <- cbind(Fst_data, q)
Fst_data_mean <- Fst_data_pop[,-c(3,1)]
Fst_data_weighted <- Fst_data_pop[,-c(3,2,1)]

# remove SVN and SVO

Fst_data_weighted <- Fst_data_weighted[-which(Fst_data_weighted$Pop1=="SVN"|Fst_data_weighted$Pop1=="SVO"|Fst_data_weighted$Pop2=="SVN"|Fst_data_weighted$Pop2=="SVO"),]

Fst_data_weighted$WeightedFst[which(Fst_data_weighted$WeightedFst>1)] <- 1
Fst_data_weighted$WeightedFst[which(Fst_data_weighted$WeightedFst<0)] <- 0

wide_weighted <- pivot_wider(Fst_data_weighted, names_from = Pop1, values_from = WeightedFst)
wide_weighted <- as.data.frame(as.matrix(wide_weighted))
rownames(wide_weighted) <- wide_weighted$Pop2
wide_weighted1 <- wide_weighted[,-c(1,2)]
wide_weighted2 <- sapply(wide_weighted1, as.numeric )
rownames(wide_weighted2) <- rownames(wide_weighted1)

#plot generally no order
# jpeg("./Figures_data/Plots/Fst_weighted_no_order.jpg", width = 2000, height = 1900)
# imagenan(wide_weighted2, lasval=2, cex.axis=6, lnumr=39, lnumc=38)
# dev.off()

# sort order by BARD
ord <- order(wide_weighted1[,6], na.last = TRUE)
wide_weighted_orderedvalues <- wide_weighted1[ord, ord]
wide_weighted_orderedvalues1 <- sapply(wide_weighted_orderedvalues, as.numeric )
rownames(wide_weighted_orderedvalues1) <- rownames(wide_weighted_orderedvalues)

#plot ordered by BARD
jpeg("./Figures_data/Plots/Fst_weighted_ordered_by_BARD.jpg", width = 3000, height = 2900)
par(mar=c(20,20,4,4))
imagenan(wide_weighted_orderedvalues1,  zlim=c(0,0.3), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6,
         col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

# order populations by long 
wide_weighted$POP <- wide_weighted$Pop2
colnames(Admix_group_dist) <-  c("POP", "Lat","Long")
weightFst_latlong <- left_join(wide_weighted, Admix_group_dist, by="POP")

rownames(weightFst_latlong) <- rownames(wide_weighted)
weightFst_latlong$Long <- as.numeric(weightFst_latlong$Long)
ord_Long <- order(weightFst_latlong$Long, na.last = TRUE)
long_sort <- wide_weighted1[ord_Long, ord_Long]
long_sort2 <- sapply(long_sort, function(x) { as.numeric(as.character(x))} )
rownames(long_sort2) <- rownames(long_sort)

#Plot ordered by Long
jpeg("./Figures_data/Plots/Fst_weighted_longsort.jpg", width = 3000, height = 2900)
par(mar=c(20,20,4,4))
imagenan(long_sort2, zlim=c(0,0.3), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6,
         col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

#--------------------------
colnames(group_distances)[1] <- "POP"
weightFst_latlong_dist <- left_join(weightFst_latlong, group_distances, by="POP")
rownames(weightFst_latlong_dist) <- rownames(wide_weighted)
weightFst_latlong_dist$BARD.y <- as.numeric(weightFst_latlong_dist$BARD.y)
ord_dist <- order(weightFst_latlong_dist$BARD.y, na.last = TRUE)
dist_sort <- wide_weighted1[ord_dist, ord_dist]
dist_sort2 <- sapply(dist_sort, function(x) { as.numeric(as.character(x))} )
rownames(dist_sort2) <- rownames(dist_sort)

#Plot ordered by dist from BARD
jpeg("./Figures_data/Plots/Fst_weighted_DistBARDsort.jpg", width = 3000, height = 2900)
par(mar=c(20,20,4,4))
imagenan(dist_sort2, zlim=c(0,0.3), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6,
         col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

#--------------------------
colnames(group_distances)[1] <- "POP"
weightFst_latlong_dist <- left_join(weightFst_latlong, group_distances, by="POP")
rownames(weightFst_latlong_dist) <- rownames(wide_weighted)
weightFst_latlong_dist$AlexNew.y <- as.numeric(weightFst_latlong_dist$AlexNew.y)
ord_dist <- order(weightFst_latlong_dist$AlexNew.y, na.last = TRUE)
dist_sort <- wide_weighted1[ord_dist, ord_dist]
dist_sort2 <- sapply(dist_sort, function(x) { as.numeric(as.character(x))} )
rownames(dist_sort2) <- rownames(dist_sort)

#Plot ordered by dist from BARD
jpeg("./Figures_data/Plots/Fst_weighted_DistAlexsort.jpg", width = 3000, height = 2900)
par(mar=c(20,20,4,4))
imagenan(dist_sort2, zlim=c(0,0.3), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6,
         col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

###################
#Fst Admixture plots

# 5 groups
# V5 - Alaska/NWT
# V2 - Europe
# V4 - Greenland
# V1 - Russia
# V3 - Saximontana

colnames(Fst_Groups) <- c("Pop", "Russia", "Europe", "Saximontana", "Greenland", "Alaska")
rownames(Fst_Groups) <- c("Russia", "Europe", "Saximontana", "Greenland", "Alaska")

#order
Fst_Groups1 <-  Fst_Groups[,-1]

ord_Fst <- order(Fst_Groups1$Alaska, na.last = TRUE)

Fst_Groups2<- Fst_Groups1[ord_Fst, ord_Fst]

# 5 groups plot
jpeg("./Figures_data/Plots/Fst_5groups.jpg", width = 2500, height = 2500)
par(mar=c(100,100,40,4))
imagenan(Fst_Groups2, zlim=c(0,0.5), yline=1,yma=15,xline=3,xma=20,lasval=2,cex.axis=6,
         col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()


########################
# Mantel test
# adjust distances and Fst to correct dimensions 

colnames(group_distances)[1] <- "POP"
weightFst_dist <- left_join(wide_weighted, group_distances, by="POP")

FST <- weightFst_dist[,3:38]
rownames(FST) <- weightFst_dist$Pop2
colnames(FST) <- weightFst_dist$Pop2

FST1<- as.data.frame(FST)
FST1$P1 <- rownames(FST)

FST2 <- FST1[,-37]
Fst_Groups_dist <- as.dist(FST2)

#---
dist <- weightFst_dist[,39:81]
rownames(dist) <- dist$POP
# compare - which don't match
dist2 <- dist %>% select(contains(".y"))
colnames(dist2) <- rownames(dist2)

dist2$P1 <- rownames(dist2)

dist3 <- dist2[,-37]
group_distances_dist <- as.dist(dist3)

#-----------------------
# not correct measure
#Fst_Groups_dist <- dist(FST2)
#group_distances_dist <- dist(dist3)

#Fst_Groups_dist[which(Fst_Groups_dist>1)] <-0.9999999
#Fst_Groups_dist[which(Fst_Groups_dist<0)] <-0

#-------------
# Mantel test
# https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
library(ade4)
mantel.rtest(group_distances_dist, Fst_Groups_dist, nrepet = 9999)

# Monte-Carlo test
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
# 
# Observation: 0.3901366 
# 
# Based on 9999 replicates
# Simulated p-value: 1e-04 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# 4.0911953253 -0.0003413571  0.0091094597


dist_Fst <- as.data.frame(cbind(Fst_Groups_dist, group_distances_dist))

jpeg("./Figures_data/Plots/Fst_Distance.jpg", width = 1000, height = 707)
plot(Fst_Groups_dist~group_distances_dist, xlab="Geographic Distance", ylab="Fst", pch=19)
dev.off()

jpeg("./Figures_data/Plots/Mantel_test.jpg", width = 3000, height = 1700)
par(mar=c(20,20,4,4))
plot((dist_Fst$Fst_Groups_dist/(1-dist_Fst$Fst_Groups_dist))~ log(dist_Fst$group_distances_dist), xlab="ln(Geographic Distance)", ylab="Fst/1-Fst", pch=19, cex = 3, mgp=c(10,5,0), cex.lab=5, cex.axis=5)
dev.off()

y <- dist_Fst$Fst_Groups_dist/(1-dist_Fst$Fst_Groups_dist)
x <- log(dist_Fst$group_distances_dist)

IBD_test <- as.data.frame(cbind(x, y))

IBD_test$x[which(is.infinite(IBD_test$x))] <- NA
model <- lm(y ~ x, IBD_test)
summary(model)

# Call:
#   lm(formula = y ~ x, data = IBD_test)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.3273 -0.2308 -0.1537 -0.0610  3.7421 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.63189    0.18036  -3.503 0.000492 ***
#   x            0.12979    0.02378   5.457 6.97e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5945 on 628 degrees of freedom
# Multiple R-squared:  0.04528,	Adjusted R-squared:  0.04376 
# F-statistic: 29.78 on 1 and 628 DF,  p-value: 6.966e-08

##############################
# plot again removing Gentian

Fst_val <- gather(FST1, key = "P2", value = "FST", 1:36)

dist_val <- gather(dist2, key = "P2", value = "Distance", 1:36)

jpeg("./Figures_data/Plots/Fst_Distance_col.jpg", width = 2000, height = 1500)
par(mar=c(10,10,4,4))
plot(Fst_val$FST~dist_val$Distance, xlab="Geographic Distance", ylab="Fst", pch=19, mgp=c(6,2,0), cex.lab=3, cex.axis=3, cex=3, col=as.factor(Fst_val$P1))
dev.off()


