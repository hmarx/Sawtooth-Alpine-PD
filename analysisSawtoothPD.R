#################################################################################################################
######## Code for: Evolutionary relationships illuminate alpine floristic diversity patterns in a remote North American wilderness
######## Load required packages and datasets ####################################################################
######## Hannah E. Marx, 18 Nov 2016 ############################################################################
#################################################################################################################

require(tidyr)
require(dplyr)
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

################################################## Taxonomy Lookup Table ################################################## 

taxonomy.table <- read.csv("data/fleshed_genera.csv")

################################################## Community Matrices ################################################## 

###### Community matrix of all plants collected ###### 
sawtooth.com.collect <- read.csv(file="data/speciespersummit_collect.csv", row.names=1)
head(sawtooth.com.collect)
split2 <- strsplit(as.character(sawtooth.com.collect$DNA.Log.Species.Name), split=" ", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.com.collect.tax <- cbind(sawtooth.com.collect[1], "genus"=genus.name2, sawtooth.com.collect[2:ncol(sawtooth.com.collect)])
head(sawtooth.com.collect.tax)
dim(sawtooth.com.collect.tax) #463 species total collected 

###### Community matrix of sequence data from MiSeq (each individual DNA accession numbers for each collection) ###### 
## MiSeq All alpine (tallus + meadow)
sawtooth.com.miseq.info <- read.csv(file="data/speciespersummit_miseqINFO.csv", row.names=1)
colnames(sawtooth.com.miseq.info) <- gsub(colnames(sawtooth.com.miseq.info), pattern = "_", replacement = " ")
sawtooth.com.miseq <- read.csv(file="data/speciespersummit_miseq.csv", row.names=2)
head(sawtooth.com.miseq)
colnames(sawtooth.com.miseq) <- gsub(colnames(sawtooth.com.miseq), pattern = "_", replacement = " ")
sawtooth.com.miseq <- sawtooth.com.miseq[2:ncol(sawtooth.com.miseq)]

## MiSeq Alpine
sawtooth.com.miseq.alpine <- read.csv(file="data/speciespersummit_miseqAlpine.csv", row.names=2)
colnames(sawtooth.com.miseq.alpine) <- gsub(colnames(sawtooth.com.miseq.alpine), pattern = "_", replacement = " ")
head(sawtooth.com.miseq.alpine)
sawtooth.com.miseq.alpine <- sawtooth.com.miseq.alpine[2:ncol(sawtooth.com.miseq.alpine)]
split2 <- strsplit(as.character(rownames(sawtooth.com.miseq.alpine)), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.com.miseq.alpine.tax <- cbind("genus"=genus.name2, sawtooth.com.miseq.alpine)
head(sawtooth.com.miseq.alpine.tax)
dim(sawtooth.com.miseq.alpine.tax) #307 alpine speices with MiSeq data

## MiSeq Meadow
sawtooth.com.miseq.meadow <- read.csv(file="data/speciespersummit_miseqMeadow.csv", row.names=2)
head(sawtooth.com.miseq.meadow)
sawtooth.com.miseq.meadow <- sawtooth.com.miseq.meadow[2:ncol(sawtooth.com.miseq.meadow)]
split2 <- strsplit(as.character(rownames(sawtooth.com.miseq.meadow)), split="_", fixed=TRUE) #split names
genus.name2 <- sapply(split2, "[", 1L)
sawtooth.com.miseq.meadow.tax <- cbind("genus"=genus.name2, sawtooth.com.miseq.meadow)
head(sawtooth.com.miseq.meadow.tax)
dim(sawtooth.com.miseq.meadow.tax) #117 meadow speices with MiSeq data ; 6 summits with meadows

###### Community Matrix Miseq for Beta phylo analyses with annotated species names and duplicate collections for each species removed (to the species level) ###### 
## MiSeq Beta All alpine (tallus + meadow)
sawtooth.com.miseq.beta <- read.csv(file="data/speciespersummit_beta.csv", row.names=2)
head(sawtooth.com.miseq.beta)
sawtooth.com.miseq.beta <- sawtooth.com.miseq.beta[2:ncol(sawtooth.com.miseq.beta)]
colnames(sawtooth.com.miseq.beta) <- gsub(colnames(sawtooth.com.miseq.beta) , pattern = "_", replacement = " ")
dim(sawtooth.com.miseq.beta) #143  unique species collected 

## MiSeq Beta Talus
sawtooth.com.miseq.talus.beta <- read.csv(file="data/speciespersummit_beta_talus.csv", row.names=2)
head(sawtooth.com.miseq.talus.beta)
sawtooth.com.miseq.talus.beta <- sawtooth.com.miseq.talus.beta[2:ncol(sawtooth.com.miseq.talus.beta)]
colnames(sawtooth.com.miseq.talus.beta) <- gsub(colnames(sawtooth.com.miseq.talus.beta), pattern = "_", replacement = " ")
dim(sawtooth.com.miseq.talus.beta) #114  unique species collected 

## MiSeq Beta Meadow
sawtooth.com.miseq.meadow.beta <- read.csv(file="data/speciespersummit_beta_meadow.csv", row.names=2)
head(sawtooth.com.miseq.meadow.beta)
sawtooth.com.miseq.meadow.beta <- sawtooth.com.miseq.meadow.beta[2:ncol(sawtooth.com.miseq.meadow.beta)]
colnames(sawtooth.com.miseq.meadow.beta) <- gsub(colnames(sawtooth.com.miseq.meadow.beta), pattern = "_", replacement = " ")
dim(sawtooth.com.miseq.meadow.beta) #74  unique species collected 


###### Community Matrix PHLAWD (to the species level) ###### 
## PHLAWD Beta All alpine (tallus + meadow)
sawtooth.com.phlawd <- read.csv(file="data/speciespersummit_phlawd.csv", row.names=2, as.is=T)
colnames(sawtooth.com.phlawd) <- gsub(colnames(sawtooth.com.phlawd) , pattern = "_", replacement = " ")
sawtooth.com.phlawd <- sawtooth.com.phlawd[2:ncol(sawtooth.com.phlawd)]
head(sawtooth.com.phlawd)
dim(sawtooth.com.phlawd) #104 alpine species with PHLAWD
split <- strsplit(as.character(rownames(sawtooth.com.phlawd)), split="_", fixed=TRUE) #split names
genus.name <- sapply(split, "[", 1L) 
sawtooth.com.phlawd.tax <- cbind("genus"=genus.name, sawtooth.com.phlawd)

## PHLAWD Beta Talus
sawtooth.com.phlawd.alpine <- read.csv(file="data/speciespersummit_phlawdAlpine.csv", row.names=2, as.is=T)
colnames(sawtooth.com.phlawd.alpine) <- gsub(colnames(sawtooth.com.phlawd.alpine) , pattern = "_", replacement = " ")
sawtooth.com.phlawd.alpine <- sawtooth.com.phlawd.alpine[2:ncol(sawtooth.com.phlawd.alpine)]
head(sawtooth.com.phlawd.alpine)
dim(sawtooth.com.phlawd.alpine) #84 alpine species with PHLAWD
split <- strsplit(as.character(rownames(sawtooth.com.phlawd.alpine)), split="_", fixed=TRUE) #split names
genus.name <- sapply(split, "[", 1L) 
sawtooth.com.phlawd.alpine.tax <- cbind("genus"=genus.name, sawtooth.com.phlawd.alpine)

## PHLAWD Beta Meadow
sawtooth.com.phlawd.meadow <- read.csv(file="data/speciespersummit_phlawdMeadow.csv", row.names=2, as.is=T)
head(sawtooth.com.phlawd.meadow)
sawtooth.com.phlawd.meadow <- sawtooth.com.phlawd.meadow[2:ncol(sawtooth.com.phlawd.meadow)]
dim(sawtooth.com.phlawd.meadow) #57 meadow species with PHLAWD; 6 summits with meadow
split <- strsplit(as.character(rownames(sawtooth.com.phlawd.meadow)), split="_", fixed=TRUE) #split names
genus.name <- sapply(split, "[", 1L)
sawtooth.com.phlawd.meadow.tax <- cbind("genus"=genus.name, sawtooth.com.phlawd.meadow)


################################################## Trees ################################################## 
########## MiSeq ML Phylogeny (individuals)
SawtoothMiseqDatedTree <- read.tree(file="output/06_Scaling/MiSeq/Sawtooth.160511.dated.bootstrap.tre") #424, dated, rooted
#write.csv(saw.phy$tip.label, "output/05_Trees/concat/MiSeq/160511/saw.phy.tips.csv")

########## MiSeq Beta Phylogeny (species)
SawtoothMiseqBeta <- read.tree("data/sawtooth.miseq.beta.tre")
SawtoothMiseqBeta$tip.label

########## PHLAWD ML Phylogeny (species)
saw.phlawd.phy <- read.nexus("output/05_Trees/concat/PHLAWD/Cipres_Data_ saw.phlawd.concat.MRE/RAxML_bipartitions.saw.phlawd.concat.MRE.nex") #105
#write.csv(saw.phlawd.phy$tip.label, "output/05_Trees/concat/PHLAWD/Cipres_Data_ saw.phlawd.concat.MRE/saw.phlawd.phy.tips.csv")


################################################## Metadata ################################################## 

# | range | elevation | WGS Lat | WGS Long
sawMeta <- read.csv("data/RAW/SawtoothMeta.csv", row.names=1)


################################################## Combine Datasets ################################################## 
dim(sawtooth.com.miseq) #424
sawMiseq <- comparative.comm(phy = SawtoothMiseqDatedTree, comm = as.matrix(t(sawtooth.com.miseq)), env = sawMeta) #as.matrix(t(sawtooth.com.miseq[, c(1,15:23)]))
sawMiseq$comm

dim(sawtooth.com.miseq.alpine) #307
sawMiseqAlp <- comparative.comm(phy = SawtoothMiseqDatedTree, comm = as.matrix(t(sawtooth.com.miseq.alpine)), env = sawMeta) #as.matrix(t(sawtooth.com.miseq[, c(1,15:23)]))
sawMiseqAlp$comm

dim(sawtooth.com.miseq.beta) #143
sawMiseqBeta <- comparative.comm(phy = SawtoothMiseqBeta, comm = as.matrix(t(sawtooth.com.miseq.beta)), env = sawMeta) #as.matrix(t(sawtooth.com.miseq[, c(1,15:23)]))
sawMiseqBeta$comm

dim(sawtooth.com.miseq.talus.beta) #114
sawMiseqBetaAlp <- comparative.comm(phy = SawtoothMiseqBeta, comm = as.matrix(t(sawtooth.com.miseq.talus.beta)), env = sawMeta) #as.matrix(t(sawtooth.com.miseq[, c(1,15:23)]))
sawMiseqBetaAlp$comm

