#################################################################################################################
######## Code for: Increasing phylogenetic stochasticity at high elevations on summits across a remote North American wilderness
######## Load required packages and datasets ####################################################################
######## Hannah E. Marx, 20 Jan 2019 ############################################################################
#################################################################################################################

require(tidyr)
require(dplyr)
require(plyr)
require(geiger)
require(ggplot2)
require(pez)
require(ape)
require(picante)
require(spacodiR)
require(pheatmap)
require(RColorBrewer)
require(phytools)
library(maps)
library(mapdata)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
require(reshape)
library(picante)
require(ecodist)
require(psych)
require(lme4)
require(MuMIn)

#rm(list = ls())
################################################## Taxonomy Lookup Table ################################################## 

taxonomy.table <- read.csv("data/fleshed_genera.csv")

################################################## Community Matrices ################################################## 

###### Community matrix of all plants collected ###### 
## Total Data 
sawtooth.com.collect <- read.csv(file="data/speciespersummit_collect.csv", row.names=1, as.is =T, check.names=FALSE, header = T)
head(sawtooth.com.collect)
dim(sawtooth.com.collect)
rownames(sawtooth.com.collect) <- sawtooth.com.collect[,1]
split2 <- strsplit(as.character(sawtooth.com.collect$Annotated_Name), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.com.collect.tax <- cbind(sawtooth.com.collect[1], "genus"=genus.name2, sawtooth.com.collect[2:ncol(sawtooth.com.collect)])
head(sawtooth.com.collect.tax)
length(unique(sawtooth.com.collect.tax$Annotated_Name)) #162 species total collected 
nrow(filter(sawtooth.com.collect.tax, Meadow == 0) %>% distinct(Annotated_Name)) #86 talus species 
nrow(filter(sawtooth.com.collect.tax, Meadow == 1) %>% distinct(Annotated_Name)) #76 Meadow species 

## Total Data Alpine
sawtooth.total.alpine <- read.csv(file="data/speciespersummit_collectAlpine.csv", row.names=2, check.names=FALSE, header = T)
head(sawtooth.total.alpine)
sawtooth.total.alpine <- sawtooth.total.alpine[2:ncol(sawtooth.total.alpine)]
split2 <- strsplit(as.character(rownames(sawtooth.total.alpine)), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.total.alpine.tax <- cbind("genus"=genus.name2, sawtooth.total.alpine)
head(sawtooth.total.alpine.tax)
dim(sawtooth.total.alpine.tax) #131 alpine speices

## Total Data Meadow
sawtooth.total.meadow <- read.csv(file="data/speciespersummit_collectMeadow.csv", row.names=2, check.names=FALSE, header = T)
head(sawtooth.total.meadow)
sawtooth.total.meadow <- sawtooth.total.meadow[2:ncol(sawtooth.total.meadow)]
split2 <- strsplit(as.character(rownames(sawtooth.total.meadow)), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.total.meadow.tax <- cbind("genus"=genus.name2, sawtooth.total.meadow)
head(sawtooth.total.meadow.tax)
dim(sawtooth.total.meadow.tax) #76 meadow speices 



################################################## Trees ################################################## 

########## Total Data ML Phylogeny (species)
saw.total.phy <- read.tree("output/09_Scaling/sawtooth.totalData.tre") 
saw.total.phy$tip.label #155


################################################## Metadata ################################################## 

# | range | elevation | WGS Lat | WGS Long
sawMeta <- read.csv("data/RAW/SawtoothMeta.csv", row.names=1)


################################################## Combine Datasets ################################################## 

#### Total Data 
dim(sawtooth.com.collect) #162
sawTotal <- comparative.comm(phy = saw.total.phy, comm = as.matrix(t(sawtooth.com.collect[,-1])), env = sawMeta) 
dim(sawTotal$comm) #154

sawTotal.alpine <- comparative.comm(phy = saw.total.phy, comm = as.matrix(t(sawtooth.total.alpine[,-1])), env = sawMeta) 
dim(sawTotal.alpine$comm) #125

sawTotal.meadow <- comparative.comm(phy = saw.total.phy, comm = as.matrix(t(sawtooth.total.meadow[,-1])), env = sawMeta) 
dim(sawTotal.meadow$comm) #73

rownames(sawtooth.com.collect)[which(!rownames(sawtooth.com.collect) %in% (saw.total.phy$tip.label))]

