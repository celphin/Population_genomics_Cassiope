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
#library(PlotSvalbard)
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
library(ggplot2)

library(SNPRelate)
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
# V5 - Greenland/NWT
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

#################################
# Mantel test - try for each group
# subset weightFst_dist by POP admix groups

# first add colnames as a row
wide_weighted[nrow(wide_weighted)+1,] <- colnames(wide_weighted)
group_distances[nrow(group_distances)+1,] <- colnames(group_distances)

# remove extra columns
wide_weighted1 <- wide_weighted[,-c(1,2)]

# determine dominant Admixture group for each population
Admix_groups <- gather(All_pop_data, Group, Amount, V1:V5, factor_key=TRUE)

Admix_groups$Group <- as.character(Admix_groups$Group)
Admix_groups <- as.data.frame(cbind(Admix_groups$Pop, Admix_groups$Group, Admix_groups$Amount))
colnames(Admix_groups) <- c("POP", "Group", "Amount")
Admix_groups$Pop <- as.factor(Admix_groups$POP)
Admix_groups$Amount <- as.numeric(Admix_groups$Amount)

maxAdmixgroup <- Admix_groups  %>%
  group_by(POP) %>%
  dplyr::summarise(MaxGroup = Group[which(Amount==max(Amount))],
                   MaxAmount=max(Amount))

maxAdmixgroups <- Admix_groups  %>%
  group_by(POP) %>%
  filter(Amount == max(Amount, na.rm=TRUE))

wide_weighted_Groups<- left_join(maxAdmixgroup, wide_weighted1, by="POP")
distances_Groups <- left_join(maxAdmixgroup, group_distances, by="POP")

#----------------------------
# split rows by admixture group for distance
distances_Groups$Group <-as.factor(distances_Groups$MaxGroup)
rows_split_dist <-  split(distances_Groups, distances_Groups$Group)

Russia_dist <- as.data.frame(rows_split_dist$V1)
# transpose all and split rows
Russia_dist[nrow(Russia_dist)+1,] <- colnames(Russia_dist)
Russia_dist <- t(as.matrix(Russia_dist))
colnames(Russia_dist) <- Russia_dist[1,]
Russia_dist <- as.data.frame(Russia_dist[-c(1:3, nrow(Russia_dist)),])
Russia_dist1 <- left_join(maxAdmixgroup, Russia_dist, by="POP")
Russia_dist2 <- Russia_dist1[which(Russia_dist1$MaxGroup=="V1"),]
Russia_dist3 <- as.data.frame(Russia_dist2[,-c(1:3)])
rownames(Russia_dist3) <- Russia_dist2$POP

Europe_dist <- as.data.frame(rows_split_dist$V2)
# transpose all and split rows
Europe_dist[nrow(Europe_dist)+1,] <- colnames(Europe_dist)
Europe_dist <- t(as.matrix(Europe_dist))
colnames(Europe_dist) <- Europe_dist[1,]
Europe_dist <- as.data.frame(Europe_dist[-c(1:3, nrow(Europe_dist)),])
Europe_dist1 <- left_join(maxAdmixgroup, Europe_dist, by="POP")
Europe_dist2 <- Europe_dist1[which(Europe_dist1$MaxGroup=="V2"),]
Europe_dist3 <- as.data.frame(Europe_dist2[,-c(1:3)])
rownames(Europe_dist3) <- Europe_dist2$POP

BC_dist <- as.data.frame(rows_split_dist$V3)
# transpose all and split rows
BC_dist[nrow(BC_dist)+1,] <- colnames(BC_dist)
BC_dist <- t(as.matrix(BC_dist))
colnames(BC_dist) <- BC_dist[1,]
BC_dist <- as.data.frame(BC_dist[-c(1:3, nrow(BC_dist)),])
BC_dist1 <- left_join(maxAdmixgroup, BC_dist, by="POP")
BC_dist2 <- BC_dist1[which(BC_dist1$MaxGroup=="V3"),]
BC_dist3 <- as.data.frame(BC_dist2[,-c(1:3)])
rownames(BC_dist3) <- BC_dist2$POP

Greenland_dist <- as.data.frame(rows_split_dist$V4)
# transpose all and split rows
Greenland_dist[nrow(Greenland_dist)+1,] <- colnames(Greenland_dist)
Greenland_dist <- t(as.matrix(Greenland_dist))
colnames(Greenland_dist) <- Greenland_dist[1,]
Greenland_dist <- as.data.frame(Greenland_dist[-c(1:3, nrow(Greenland_dist)),])
Greenland_dist1 <- left_join(maxAdmixgroup, Greenland_dist, by="POP")
Greenland_dist2 <- Greenland_dist1[which(Greenland_dist1$MaxGroup=="V4"),]
Greenland_dist3 <- as.data.frame(Greenland_dist2[,-c(1:3)])
rownames(Greenland_dist3) <- Greenland_dist2$POP

Alaska_dist <- as.data.frame(rows_split_dist$V5)
# transpose all and split rows
Alaska_dist[nrow(Alaska_dist)+1,] <- colnames(Alaska_dist)
Alaska_dist <- t(as.matrix(Alaska_dist))
colnames(Alaska_dist) <- Alaska_dist[1,]
Alaska_dist <- as.data.frame(Alaska_dist[-c(1:3, nrow(Alaska_dist)),])
Alaska_dist1 <- left_join(maxAdmixgroup, Alaska_dist, by="POP")
Alaska_dist2 <- Alaska_dist1[which(Alaska_dist1$MaxGroup=="V5"),]
Alaska_dist3 <- as.data.frame(Alaska_dist2[,-c(1:3)])
rownames(Alaska_dist3) <- Alaska_dist2$POP

#----------------------------
# split rows by admixture group for distance
wide_weighted_Groups$Group <-as.factor(wide_weighted_Groups$MaxGroup)
rows_split_Fst <-  split(wide_weighted_Groups, wide_weighted_Groups$Group)


Russia_Fst <- as.data.frame(rows_split_Fst$V1)
# transpose all and split rows
Russia_Fst[nrow(Russia_Fst)+1,] <- colnames(Russia_Fst)
Russia_Fst <- t(as.matrix(Russia_Fst))
colnames(Russia_Fst) <- Russia_Fst[1,]
Russia_Fst <- as.data.frame(Russia_Fst[-c(1:3, nrow(Russia_Fst)),])
Russia_Fst1 <- left_join(maxAdmixgroup, Russia_Fst, by="POP")
Russia_Fst2 <- Russia_Fst1[which(Russia_Fst1$MaxGroup=="V1"),]
Russia_Fst3 <- as.data.frame(Russia_Fst2[,-c(1:3)])
rownames(Russia_Fst3) <- Russia_Fst2$POP

Europe_Fst <- as.data.frame(rows_split_Fst$V2)
# transpose all and split rows
Europe_Fst[nrow(Europe_Fst)+1,] <- colnames(Europe_Fst)
Europe_Fst <- t(as.matrix(Europe_Fst))
colnames(Europe_Fst) <- Europe_Fst[1,]
Europe_Fst <- as.data.frame(Europe_Fst[-c(1:3, nrow(Europe_Fst)),])
Europe_Fst1 <- left_join(maxAdmixgroup, Europe_Fst, by="POP")
Europe_Fst2 <- Europe_Fst1[which(Europe_Fst1$MaxGroup=="V2"),]
Europe_Fst3 <- as.data.frame(Europe_Fst2[,-c(1:3)])
rownames(Europe_Fst3) <- Europe_Fst2$POP

BC_Fst <- as.data.frame(rows_split_Fst$V3)
# transpose all and split rows
BC_Fst[nrow(BC_Fst)+1,] <- colnames(BC_Fst)
BC_Fst <- t(as.matrix(BC_Fst))
colnames(BC_Fst) <- BC_Fst[1,]
BC_Fst <- as.data.frame(BC_Fst[-c(1:3, nrow(BC_Fst)),])
BC_Fst1 <- left_join(maxAdmixgroup, BC_Fst, by="POP")
BC_Fst2 <- BC_Fst1[which(BC_Fst1$MaxGroup=="V3"),]
BC_Fst3 <- as.data.frame(BC_Fst2[,-c(1:3)])
rownames(BC_Fst3) <- BC_Fst2$POP

Greenland_Fst <- as.data.frame(rows_split_Fst$V4)
# transpose all and split rows
Greenland_Fst[nrow(Greenland_Fst)+1,] <- colnames(Greenland_Fst)
Greenland_Fst <- t(as.matrix(Greenland_Fst))
colnames(Greenland_Fst) <- Greenland_Fst[1,]
Greenland_Fst <- as.data.frame(Greenland_Fst[-c(1:3, nrow(Greenland_Fst)),])
Greenland_Fst1 <- left_join(maxAdmixgroup, Greenland_Fst, by="POP")
Greenland_Fst2 <- Greenland_Fst1[which(Greenland_Fst1$MaxGroup=="V4"),]
Greenland_Fst3 <- as.data.frame(Greenland_Fst2[,-c(1:3)])
rownames(Greenland_Fst3) <- Greenland_Fst2$POP

Alaska_Fst <- as.data.frame(rows_split_Fst$V5)
# transpose all and split rows
Alaska_Fst[nrow(Alaska_Fst)+1,] <- colnames(Alaska_Fst)
Alaska_Fst <- t(as.matrix(Alaska_Fst))
colnames(Alaska_Fst) <- Alaska_Fst[1,]
Alaska_Fst <- as.data.frame(Alaska_Fst[-c(1:3, nrow(Alaska_Fst)),])
Alaska_Fst1 <- left_join(maxAdmixgroup, Alaska_Fst, by="POP")
Alaska_Fst2 <- Alaska_Fst1[which(Alaska_Fst1$MaxGroup=="V5"),]
Alaska_Fst3 <- as.data.frame(Alaska_Fst2[,-c(1:3)])
rownames(Alaska_Fst3) <- Alaska_Fst2$POP

#--------------------------
# # rerun mantel test for all 5 groups

# https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
library(ade4)

#---------------------
# Alaska
Alaska_Fst3 <- Alaska_Fst3[-5,-5]
Alaska_dist3 <- Alaska_dist3[-5,-5]

Pop <- as.data.frame(expand_grid(rownames(Alaska_Fst3), rownames(Alaska_Fst3)))
Pop_nam <- paste0(Pop$`rownames(Alaska_Fst3)...1`, "_", Pop$`rownames(Alaska_Fst3)...2`)
  
Alaska_Fst4 <- as.dist(as.matrix(Alaska_Fst3))
Alaska_dist4 <- as.dist(as.matrix(Alaska_dist3))
mantel.rtest(Alaska_dist4, Alaska_Fst4, nrepet = 9999)

summary(lm(Alaska_Fst4~Alaska_dist4))

#Observation: 0.6421697 
#Simulated p-value: 1e-04 
#Std.Obs  Expectation     Variance 
#3.564626859 -0.001881282  0.032644674 

png("./Figures_data/Plots/Fst_Distance_Alaska.png", width = 1000, height = 707)
par(mar=c(7,7,4,4))
plot(Alaska_Fst4~Alaska_dist4,  xlab="Geographic Distance (km)", ylab="Fst", pch=19, cex = 1.5, cex.lab=2, cex.axis=2)
dev.off()

png("./Figures_data/Plots/Mantel_test_Alaska.jpg", width = 3000, height = 1700)
par(mar=c(20,20,4,4))
plot((Alaska_Fst4/(1-Alaska_Fst4))~ log(Alaska_dist4), xlab="ln(Geographic Distance)", ylab="Fst/1-Fst", pch=19, cex = 3, mgp=c(10,5,0), cex.lab=5, cex.axis=5)
dev.off()

# add labels to plot
Pop_nam1 <- Pop_nam[1:length(Pop_nam)/2]

jpeg("./Figures_data/Plots/Fst_Distance_Alaska_labels.jpg", width = 1000, height = 707)
plot(Alaska_Fst4[1:length(Pop_nam)/2]~Alaska_dist4[1:length(Pop_nam)/2],  xlab="Geographic Distance", ylab="Fst", pch=19)
text(Alaska_Fst4[1:length(Pop_nam)/2]~Alaska_dist4[1:length(Pop_nam)/2], labels=Pop_nam1, cex=0.9, font=2)
dev.off()

y <- (Alaska_Fst4/(1-Alaska_Fst4))
x <- log(Alaska_dist4)

IBD_test <- as.data.frame(cbind(x, y))

IBD_test$x[which(is.infinite(IBD_test$x))] <- NA
model <- lm(y ~ x, IBD_test)
summary(model)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.08848    0.02849  -3.106  0.00238 ** 
#   x            0.03309    0.00419   7.897 1.65e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05179 on 118 degrees of freedom
# Multiple R-squared:  0.3458,	Adjusted R-squared:  0.3402 
# F-statistic: 62.37 on 1 and 118 DF,  p-value: 1.646e-12

###################################################
# Europe
Europe_Fst3 <- Europe_Fst3[-9,-9]
Europe_dist3 <- Europe_dist3[-9,-9]

Pop <- as.data.frame(expand_grid(rownames(Europe_Fst3), rownames(Europe_Fst3)))
Pop_nam <- paste0(Pop$`rownames(Europe_Fst3)...1`, "_", Pop$`rownames(Europe_Fst3)...2`)

Europe_Fst4 <- as.dist(as.matrix(Europe_Fst3))
Europe_dist4 <- as.dist(as.matrix(Europe_dist3))
mantel.rtest(Europe_dist4, Europe_Fst4, nrepet = 9999)

summary(lm(Europe_Fst4~Europe_dist4))

# Observation: 0.7593612 
# Simulated p-value: 0.0012 
# Std.Obs  Expectation     Variance 
# 5.0255499191 0.0009161117 0.0227761940  

png("./Figures_data/Plots/Fst_Distance_Europe.png", width = 1000, height = 707)
par(mar=c(7,7,4,4))
plot(Europe_Fst4~Europe_dist4,  xlab="Geographic Distance (km)", ylab="Fst", pch=19, cex = 1.5, cex.lab=2, cex.axis=2)
dev.off()

png("./Figures_data/Plots/Mantel_test_Europe.png", width = 3000, height = 1700)
par(mar=c(20,20,4,4))
plot((Europe_Fst4/(1-Europe_Fst4))~ log(Europe_dist4), xlab="ln(Geographic Distance)", ylab="Fst/1-Fst", pch=19, cex = 3, mgp=c(10,5,0), cex.lab=5, cex.axis=5)
dev.off()

# add labels to plot
Pop_nam1 <- Pop_nam[1:length(Pop_nam)/2]

jpeg("./Figures_data/Plots/Fst_Distance_Europe_labels.jpg", width = 1000, height = 707)
plot(Europe_Fst4[1:length(Pop_nam)/2]~Europe_dist4[1:length(Pop_nam)/2],  xlab="Geographic Distance", ylab="Fst", pch=19)
text(Europe_Fst4[1:length(Pop_nam)/2]~Europe_dist4[1:length(Pop_nam)/2], labels=Pop_nam1, cex=0.9, font=2)
dev.off()

y <- (Europe_Fst4/(1-Europe_Fst4))
x <- log(Europe_dist4)

IBD_test <- as.data.frame(cbind(x, y))

IBD_test$x[which(is.infinite(IBD_test$x))] <- NA
model <- lm(y ~ x, IBD_test)
summary(model)

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.075681 -0.023191  0.004281  0.017179  0.072228 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.021994   0.024627  -0.893    0.377    
# x            0.020631   0.003487   5.917 4.83e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03482 on 43 degrees of freedom
# Multiple R-squared:  0.4488,	Adjusted R-squared:  0.436 
# F-statistic: 35.01 on 1 and 43 DF,  p-value: 4.832e-07

###################################################
# Greenland
Greenland_Fst3 <- Greenland_Fst3[-5,-5]
Greenland_dist3 <- Greenland_dist3[-5,-5]

Pop <- as.data.frame(expand_grid(rownames(Greenland_Fst3), rownames(Greenland_Fst3)))
Pop_nam <- paste0(Pop$`rownames(Greenland_Fst3)...1`, "_", Pop$`rownames(Greenland_Fst3)...2`)

Greenland_Fst4 <- as.dist(as.matrix(Greenland_Fst3))
Greenland_dist4 <- as.dist(as.matrix(Greenland_dist3))
mantel.rtest(Greenland_dist4, Greenland_Fst4, nrepet = 9999)

#Observation: -0.2967934 
#Simulated p-value: 0.8697 
#Std.Obs  Expectation     Variance 
#-0.958513472  0.002183467  0.097292382

jpeg("./Figures_data/Plots/Fst_Distance_Greenland.jpg", width = 1000, height = 707)
plot(Greenland_Fst4~Greenland_dist4,  xlab="Geographic Distance", ylab="Fst", pch=19)
dev.off()

jpeg("./Figures_data/Plots/Mantel_test_Greenland.jpg", width = 3000, height = 1700)
par(mar=c(20,20,4,4))
plot((Greenland_Fst4/(1-Greenland_Fst4))~ log(Greenland_dist4), xlab="ln(Geographic Distance)", ylab="Fst/1-Fst", pch=19, cex = 3, mgp=c(10,5,0), cex.lab=5, cex.axis=5)
dev.off()

# add labels to plot
Pop_nam1 <- Pop_nam[1:length(Pop_nam)/2]

jpeg("./Figures_data/Plots/Fst_Distance_Greenland_labels.jpg", width = 1000, height = 707)
plot(Greenland_Fst4[1:length(Pop_nam)/2]~Greenland_dist4[1:length(Pop_nam)/2],  xlab="Geographic Distance", ylab="Fst", pch=19)
text(Greenland_Fst4[1:length(Pop_nam)/2]~Greenland_dist4[1:length(Pop_nam)/2], labels=Pop_nam1, cex=0.9, font=2)
dev.off()

y <- (Greenland_Fst4/(1-Greenland_Fst4))
x <- log(Greenland_dist4)

IBD_test <- as.data.frame(cbind(x, y))

IBD_test$x[which(is.infinite(IBD_test$x))] <- NA
model <- lm(y ~ x, IBD_test)
summary(model)


# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.11093 -0.06886 -0.02355  0.01532  0.18416 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.33279    0.09427   3.530  0.00224 **
#   x           -0.03109    0.01563  -1.989  0.06128 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0908 on 19 degrees of freedom
# Multiple R-squared:  0.1724,	Adjusted R-squared:  0.1288 
# F-statistic: 3.957 on 1 and 19 DF,  p-value: 0.06128

###################################################


