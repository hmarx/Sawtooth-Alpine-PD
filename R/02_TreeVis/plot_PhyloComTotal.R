#################################################################################################################
######## Plot Sawtooth Community Phylogeny   ####################################################################
######## Load required packages and datasets ####################################################################
######## Hannah E. Marx, 20 Jan 2019 ############################################################################
#################################################################################################################


taxonomyDir <- "output/09_Scaling/Congruify/fleshed_genera.csv"
refDates <- "output/09_Scaling/Congruify/out_dates.tre"

source("R/functions/treeFunctions.R")
source("R/functions/plot.colorCom.R")
source("R/functions/trait.plot.colorTips.R")


########################### Plot Tree ###########################################
###### Color taxa on summit: Code modified from picante color.plot.phylo
##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy (dated with PATHd8)

rownames(sawTotal$comm)
## Add a row of species names
com.data.plot <- as.data.frame(t(sawTotal$comm))
head(com.data.plot)

## Re-order column names (increasing elevation)
com.data.plot <- com.data.plot[,c("Horstmann Peak",
                                  "Braxon Peak", 
                                  "Thompson Peak",
                                  "Snowyside Peak",
                                  "Mount Cramer",
                                  "D.O. Lee Peak",
                                  "Salzburger Spitzl",
                                  "Castle Peak",
                                  "Hyndman Peak")]
                                  
head(com.data.plot)
# convert to numerics
com.data.plot[2:ncol(com.data.plot)] <- sapply(com.data.plot[2:ncol(com.data.plot)] , 
                                               function(x) as.numeric(as.character(x)))

cols = list("Horstmann Peak"= c("white", "#E76BF3"),
            "Braxon Peak"= c("white", "#FF67A4"),
            "Thompson Peak"= c("white", "#00C0AF"),
            "Snowyside Peak"= c("white", "#B983FF"),
            "Mount Cramer" = c("white", "#00B0F6"),
            "D.O. Lee Peak" = c("white", "#E58700"),
            "Salzburger Spitzl" = c("white", "#A3A500"),
            "Castle Peak" = c("white", "#00BA38"),
            "Hyndman Peak" = c("white", "#F8766D"))

sawTotalplot <- sawTotal$phy


####################### ALl Alpine  (Tallus + Meadow) with "Meadow" labeled
head(splist.raw.mi.tmp)

splist.meadow <- unique(na.omit(sawtooth.com.collect[sawtooth.com.collect$Meadow == 1, c(1,2)]))
splist.meadow

com.total.plot.meadow2 <- merge(splist.meadow, t(sawTotal$comm), by.x=1, by.y =0, all.x=F, all.y=T)
com.total.plot.meadow2[is.na(com.total.plot.meadow2[,2]),2] <- 0
rownames(com.total.plot.meadow2) <- com.total.plot.meadow2$Accepted.Name.GenBank
#com.total.plot.meadow2 <- com.total.plot.meadow2[rownames(com.total.plot.meadow2) != "Castilleja_applegatei", ]

com.total.plot.meadow3 <- cbind(GenBank = rep(0), MiSeq = rep(0), com.total.plot.meadow2)
head(com.total.plot.meadow3) 

accession.collections <- read.csv("figs/tables/accession.collections.csv", row.names = 1)

accession.collections.df <- gather(accession.collections[1:8], "data", "gene", 3:8)
tmpdata <- strsplit(accession.collections.df$gene, split="|", fixed=TRUE)
data.name <- sapply(tmpdata, "[", 1L) 
accession.collections.df <- na.omit(cbind(accession.collections.df[1], data.name))
accession.collections.df <- distinct(accession.collections.df)

Miseqtax <- as.character(accession.collections.df[accession.collections.df$data.name == "miseq",1])
GenBanktax <- as.character(accession.collections.df[accession.collections.df$data.name == "ncbi",1])
 
com.total.plot.meadow3$MiSeq[com.total.plot.meadow3$Accepted.Name.GenBank %in% Miseqtax] <- 1
com.total.plot.meadow3$GenBank[com.total.plot.meadow3$Accepted.Name.GenBank %in% GenBanktax] <- 1
com.total.plot.meadow3 <- com.total.plot.meadow3
head(com.total.plot.meadow3)
com.total.plot.meadow3$Accepted.Name.GenBank <- gsub(com.total.plot.meadow3$Accepted.Name.GenBank, pattern = "var_", replacement = "var._")
com.total.plot.meadow3$Accepted.Name.GenBank <- gsub(com.total.plot.meadow3$Accepted.Name.GenBank, pattern = "ssp_", replacement = "ssp._")
rownames(com.total.plot.meadow3) <- gsub(rownames(com.total.plot.meadow3), pattern = "var_", replacement = "var._")
rownames(com.total.plot.meadow3) <- gsub(rownames(com.total.plot.meadow3), pattern = "ssp_", replacement = "ssp._")

#rownames(com.total.plot.meadow3[rownames(com.total.plot.meadow3) =="Draba_crassifolia",]) <- "Draba_sp."
#rownames(com.total.plot.meadow3[rownames(com.total.plot.meadow3) =="Phlox",]) <- "Phlox_sp."
#rownames(com.total.plot.meadow3[rownames(com.total.plot.meadow3) =="Erigeron",]) <- "Erigeron_sp."
#rownames(com.total.plot.meadow3)[rownames(com.total.plot.meadow3) =="Draba_crassifolia"] <- "Draba_sp."
#rownames(com.total.plot.meadow3)[rownames(com.total.plot.meadow3) =="Phlox"] <- "Phlox_sp."
#rownames(com.total.plot.meadow3)[rownames(com.total.plot.meadow3) =="Erigeron"] <- "Erigeron_sp."

#com.total.plot.meadow3$Accepted.Name.GenBank[rownames(com.total.plot.meadow3) =="Draba_sp."] <- "Draba_sp."
#com.total.plot.meadow3$Accepted.Name.GenBank[rownames(com.total.plot.meadow3) =="Phlox_sp."] <- "Phlox_sp."
#com.total.plot.meadow3$Accepted.Name.GenBank[rownames(com.total.plot.meadow3) =="Erigeron_sp."] <- "Erigeron_sp."

sawTotalplot$tip.label <- gsub(sawTotalplot$tip.label, pattern = "var_", replacement = "var._")
sawTotalplot$tip.label <- gsub(sawTotalplot$tip.label, pattern = "ssp_", replacement = "ssp._")

#sawTotalplot$tip.label[sawTotalplot$tip.label== "Draba_crassifolia"] <- "Draba_sp."
#sawTotalplot$tip.label[sawTotalplot$tip.label== "Phlox"] <- "Phlox_sp."
#sawTotalplot$tip.label[sawTotalplot$tip.label== "Erigeron"] <- "Erigeron_sp."

pdf("figs/tree/SawtoothCommunityPhyloBetaColorMeadow.pdf", 10,10)
trait.plot.colorTip(tree = sawTotalplot,
                    dat = com.total.plot.meadow3,
                    shape1 = "GenBank",
                    pch1 = 8,
                    shape2 = "MiSeq",
                    pch2 = 18,
                    cex.shape = 0.75,
                    cols = list("Horstmann Peak"= c("white", "#E76BF3"),
                                "Braxon Peak"= c("white", "#FF67A4"),
                                "Thompson Peak"= c("white", "#00C0AF"),
                                "Snowyside Peak"= c("white", "#B983FF"),
                                "Mount Cramer" = c("white", "#00B0F6"),
                                "D.O. Lee Peak" = c("white", "#E58700"),
                                "Salzburger Spitzl" = c("white", "#A3A500"),
                                "Castle Peak" = c("white", "#00BA38"),
                                "Hyndman Peak" = c("white", "#F8766D"),
                                "Meadow" = c("white", "grey"),
                                "GenBank" = c("white", "black"),
                                "MiSeq" = c("white", "black")),
                    col.names = c("black", "black"),
                    trait = "Meadow",
                    datTr = com.total.plot.meadow3,
                    taxa.names = "Accepted.Name.GenBank",
                    cex.lab = 0.4,
                    font.lab = 1)
dev.off()

