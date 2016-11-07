

########################################################### 6_Visualize Tree ######################################################### 

source("R/treeFunctions.R")
source("R/plot.colorCom.R")
source("R/trait.plot.colorTips.R")

SawtoothMiseqDatedTree <- read.tree(file="output/06_Scaling/MiSeq/Sawtooth.160511.dated.bootstrap.tre") #424 (alpine & meadow), dated, rooted
sawMiseqBeta$phy # 143; Just alpine species


########################### Plot Tree ###########################################
###### Color taxa on summit: Code modified from picante color.plot.phylo
##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy (dated with PATHd8)

####################### ALl Alpine  (Tallus + Meadow)
######## Plot phylo with p/a on summits plotted across tips
rownames(sawMiseqBeta$comm)
## Add a row of species names
com.data.plot <- as.data.frame(t(sawMiseqBeta$comm))
head(com.data.plot)

## Re-order column names (decreasing species richness, MiSeq seqs)
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

cols = list("Hyndman Peak" = c("white", "#F8766D"),
            "D.O. Lee Peak" = c("white", "#E58700"),
            "Salzburger Spitzl" = c("white", "#A3A500"),
            "Castle Peak" = c("white", "#00BA38"),
            "Thompson Peak"= c("white", "#00C0AF"),
            "Mount Cramer" = c("white", "#00B0F6"),
            "Snowyside Peak"= c("white", "#B983FF"),
            "Horstmann Peak"= c("white", "#E76BF3"),
            "Braxon Peak"= c("white", "#FF67A4"))

# plot
pdf("output/07_VisualizeTree/MiSeq/SawtoothCommunityPhyloBetaColor.pdf", 10,10)
plot.colorCom(tree = sawMiseqBeta$phy,
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


####################### ALl Alpine  (Tallus + Meadow) with "Meadow" labeled
head(sawtooth.com.miseq.info)

splist.meadow <- unique(na.omit(splist.raw.mi.tmp[splist.raw.mi.tmp$Meadow == 1, c(1,13)]))
splist.meadow

com.data.plot.meadow2 <- merge(splist.meadow, t(sawMiseqBeta$comm), by.x=1, by.y =0, all.x=F, all.y=T)
com.data.plot.meadow2[is.na(com.data.plot.meadow2[,2]),2] <- 0
rownames(com.data.plot.meadow2) <- com.data.plot.meadow2$Annotated.Name
 
pdf("output/07_VisualizeTree/MiSeq/SawtoothCommunityPhyloTips.BetaColor.pdf", 10,10)
trait.plot.colorTip(tree = sawMiseqBeta$phy,
                    dat = com.data.plot.meadow2,
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
                    datTr = com.data.plot.meadow2,
                    taxa.names = "Annotated.Name",
                    col.names = c("black", "black"),
                    cex.lab = 0.5,
                    font.lab = 1)
dev.off()


cols = list("Castle Peak" = c("white", "red4"),
            "Hyndman Peak" = c("white", "darkorange1"),
            "Salzburger Spitzl" = c("white", "yellow"),
            "Mount Cramer" = c("white", "chartreuse"),
            "D.O. Lee Peak" = c("white", "goldenrod"),
            "Horstmann Peak"= c("white", "darkgreen"),
            "Braxon Peak"= c("white", "purple2"),
            "Snowyside Peak"= c("white", "hotpink"),
            "Thompson Peak"= c("white", "deeppink4"),
            "Meadow" = c("white", "grey")),
