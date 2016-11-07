

########################################################### 6_Visualize Tree ######################################################### 

source("R/treeFunctions.R")
source("R/plot.colorTips.R")
source("R/trait.plot.colorTips.R")

SawtoothMiseqDatedTree <- read.tree(file="output/06_Scaling/MiSeq/Sawtooth.160511.dated.bootstrap.tre") #424 (alpine & meadow), dated, rooted
sawMiseqAlp$phy # 307; Just alpine species


########################### Plot Tree ###########################################
###### Color taxa on summit: Code modified from picante color.plot.phylo
##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy (dated with PATHd8)

####################### ALl Alpine  (Tallus + Meadow)
######## Plot phylo with p/a on summits plotted across tips
rownames(sawMiseq$comm)
## Add a row of species names
com.data.plot <- as.data.frame(t(sawMiseq$comm))
head(com.data.plot)

## Re-order column names (increasing species richness)
com.data.plot <- com.data.plot[,c("Hyndman Peak", 
                                  "D.O. Lee Peak",
                                  "Salzburger Spitzl",
                                  "Castle Peak",
                                  "Thompson Peak",
                                  "Mount Cramer",
                                  "Snowyside Peak",
                                  "Horstmann Peak",
                                  "Braxon Peak")]
head(com.data.plot)
# convert to numerics
com.data.plot[2:ncol(com.data.plot)] <- sapply(com.data.plot[2:ncol(com.data.plot)] , 
                                               function(x) as.numeric(as.character(x)))

# plot
pdf("output/07_VisualizeTree/MiSeq/SawtoothCommunityPhylo.pdf", 10,10)
plot.colorCom(tree = sawMiseq$phy,
              dat = com.data.plot,
              cols = list("Hyndman Peak" = c("white", "#F8766D"),
                          "D.O. Lee Peak" = c("white", "#E58700"),
                          "Salzburger Spitzl" = c("white", "#A3A500"),
                          "Castle Peak" = c("white", "#00BA38"),
                          "Thompson Peak"= c("white", "#00C0AF"),
                          "Mount Cramer" = c("white", "#00B0F6"),
                          "Snowyside Peak"= c("white", "#B983FF"),
                          "Horstmann Peak"= c("white", "#E76BF3"),
                          "Braxon Peak"= c("white", "#FF67A4")),
              col.names = c("grey", "black"),
              cex.lab = 0.4,
              font.lab = 0.5,
              w = 1/20)

dev.off()


####################### Just Alpine Tallus
######## Plot phylo with p/a on summits plotted across tips
rownames(sawMiseqAlp$comm)
## Add a row of species names
com.data.plot.alpine <- as.data.frame(t(sawMiseqAlp$comm))
head(com.data.plot.alpine)

## Re-order column names (increasing species richness)
com.data.plot.alpine <- com.data.plot.alpine[,c("Hyndman Peak", 
                                             "D.O. Lee Peak",
                                             "Salzburger Spitzl",
                                             "Castle Peak",
                                             "Thompson Peak",
                                             "Mount Cramer",
                                             "Snowyside Peak",
                                             "Horstmann Peak",
                                             "Braxon Peak")]
head(com.data.plot.alpine)
# convert to numerics
com.data.plot.alpine[2:ncol(com.data.plot.alpine)] <- sapply(com.data.plot.alpine[2:ncol(com.data.plot.alpine)] , 
                                               function(x) as.numeric(as.character(x)))

# plot
pdf("output/07_VisualizeTree/MiSeq/SawtoothCommunityPhyloAlpine.pdf", 10,10)
plot.colorCom(tree = sawMiseqAlp$phy,
                    dat = com.data.plot.alpine,
                    cols = list("Castle Peak" = c("white", "red4"),
                                "Hyndman Peak" = c("white", "darkorange1"),
                                "Salzburger Spitzl" = c("white", "yellow"),
                                "Mount Cramer" = c("white", "chartreuse"),
                                "D.O. Lee Peak" = c("white", "goldenrod"),
                                "Horstmann Peak"= c("white", "darkgreen"),
                                "Braxon Peak"= c("white", "purple2"),
                                "Snowyside Peak"= c("white", "hotpink"),
                                "Thompson Peak"= c("white", "deeppink4")),
                    col.names = c("grey", "black"),
                    cex.lab = 0.4,
                    font.lab = 0.5,
                    w = 1/20)

dev.off()



head(sawtooth.com.miseq.info)
com.data.plot.meadow <- as.data.frame(cbind(species = sawtooth.com.miseq.info$Tip.Label.160511, sawtooth.com.miseq.info[c(10, 16:24)]))
head(com.data.plot.meadow)
com.data.plot.meadow <- com.data.plot.meadow[!is.na(com.data.plot.meadow$species),]
# convert to numerics
com.data.plot.meadow[2:ncol(com.data.plot.meadow)] <- sapply(com.data.plot.meadow[2:ncol(com.data.plot.meadow)] , 
                                               function(x) as.numeric(as.character(x)))
rownames(com.data.plot.meadow) <- com.data.plot.meadow$species
treedata(sawMiseq$phy, com.data.plot.meadow)

## Re-order column names (increasing species richness)
com.data.plot.meadow <- com.data.plot.meadow[,c("species", "Meadow", "Hyndman Peak", 
                                                "D.O. Lee Peak",
                                                "Salzburger Spitzl",
                                                "Castle Peak",
                                                "Thompson Peak",
                                                "Mount Cramer",
                                                "Snowyside Peak",
                                                "Horstmann Peak",
                                                "Braxon Peak")]
colnames(com.data.plot.meadow)


pdf("output/07_VisualizeTree/MiSeq/SawtoothCommunityPhyloTipsColor.pdf", 10,10)
trait.plot.colorTip(tree = sawMiseq$phy,
                    dat = com.data.plot.meadow,
                    cols = list("Hyndman Peak" = c("white", "#F8766D"),
                                "D.O. Lee Peak" = c("white", "#E58700"),
                                "Salzburger Spitzl" = c("white", "#A3A500"),
                                "Castle Peak" = c("white", "#00BA38"),
                                "Thompson Peak"= c("white", "#00C0AF"),
                                "Mount Cramer" = c("white", "#00B0F6"),
                                "Snowyside Peak"= c("white", "#B983FF"),
                                "Horstmann Peak"= c("white", "#E76BF3"),
                                "Braxon Peak"= c("white", "#FF67A4"),
                                "Meadow" = c("white", "grey")),
                    trait = "Meadow",
                    datTr = com.data.plot.meadow,
                    taxa.names = "species",
                    col.names = c("black", "black"),
                    cex.lab = 0.25,
                    font.lab = 0.5)
dev.off()
