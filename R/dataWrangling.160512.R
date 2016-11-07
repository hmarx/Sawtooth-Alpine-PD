#### Data wranagling 

################################################## MiSeq Processing ################################################## 

############ Correct 141219 Sample Sheet (for processesing MiSeq reads)  ############ 
collections <- read.csv("../sequences/runs/141219/SawtoothSamples_extractionCollectionInfoV2sub.csv")
sample141219 <- read.csv("../sequences/runs/141219/SubmissionTMP.csv")
head(collections)
head(sample141219)

sample141219.new <- merge(sample141219, collections, by.x=1, by.y=3, all.x=T, all.y=F)
#write.csv(sample141219.new, "../sequences/runs/141219/SubmissionTMPmodified.csv")

################################################## Community Matrix ################################################## 

splist.raw.mi.tmp <- read.csv("data/RAW/SawtoothCollectionInfo.csv", as.is=T)
head(splist.raw.mi.tmp)
dim(splist.raw.mi.tmp) #489

splist.raw.mi.tmp.alp <- splist.raw.mi.tmp[splist.raw.mi.tmp$Meadow != 1,] #true alpine species 
dim(splist.raw.mi.tmp.alp) # 349 talus species
splist.raw.mi.tmp.mead <- splist.raw.mi.tmp[splist.raw.mi.tmp$Meadow == 1,] #meadow species
dim(splist.raw.mi.tmp.mead) #140 talus species 

head(splist.raw.mi.tmp)

########################## Summarize Collections, MiSeq and PHLAWD

sum1 <- splist.raw.mi.tmp %>% group_by(Summit) %>% summarise(n_collected = n_distinct(Collection.. ), n_MiSeq = n_distinct(Tip.Label.160511),
                                                     n_PHLAWD = n_distinct(PHLAWD))

sum2 <- filter(splist.raw.mi.tmp, Meadow == 0) %>%  group_by(Summit) %>% summarise(n_talus_collected = n_distinct(Collection.. ), 
                                                                                   n_talus_MiSeq = n_distinct(Tip.Label.160511),
                                                                                   n_talus_PHLAWD = n_distinct(PHLAWD))

sum3 <- filter(splist.raw.mi.tmp, Meadow == 1) %>%  group_by(Summit) %>% summarise(n_meadow_collected = n_distinct(Collection.. ), 
                                                                           n_meadow_MiSeq = n_distinct(Tip.Label.160511),
                                                                           n_meadow_PHLAWD = n_distinct(PHLAWD))

sumtot <- merge(cbind(sum1, sum2[-1]), sum3, by =1, all.x=T)
write.csv(sumtot, file="output/summary.collection.csv")

############ All Species collected  ############ 

splist.collect.alpine.info <- splist.raw.mi.tmp[!is.na(splist.raw.mi.tmp.alp$DNA.Log.Species.Name), c(5, 10, 11)] 
head(splist.collect.alpine.info)
splist.col <- splist.collect.alpine.info %>% group_by(Summit) %>% distinct(DNA.Log.Species.Name) 

#### Talus
saw.com.collect.alpine <- spread(data = splist.col, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.collect.alpine) #138 alpine species collected 
head(saw.com.collect.alpine)
colnames(saw.com.collect.alpine) <- gsub(colnames(saw.com.collect.alpine), pattern = " ", replacement = "_")
write.csv(saw.com.collect.alpine, "data/speciespersummit_collectAlpine.csv")

#### Meadow
splist.collect.meadow.info <- splist.raw.mi.tmp.mead[!is.na(splist.raw.mi.tmp.mead$DNA.Log.Species.Name), c(5, 10, 11)] 
head(splist.collect.meadow.info)
splist.col <- splist.collect.meadow.info %>% group_by(Summit) %>% distinct(DNA.Log.Species.Name) 

saw.com.collect.meadow <- spread(data = splist.col, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.collect.meadow) #79 meadow species collected 
head(saw.com.collect.meadow)
colnames(saw.com.collect.meadow) <- gsub(colnames(saw.com.collect.meadow), pattern = " ", replacement = "_")
write.csv(saw.com.collect.meadow, "data/speciespersummit_collectMeadow.csv")

############ Community Matrix MiSeq ############ 

splist.col <- splist.raw.mi.tmp[!is.na(splist.raw.mi.tmp$Tip.Label.160511),] #424

saw.com.collect.alpine <- spread(data = splist.col, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.collect.alpine) #432 alpine species collected 
head(saw.com.collect.alpine)
colnames(saw.com.collect.alpine) <- gsub(colnames(saw.com.collect.alpine), pattern = " ", replacement = "_")
write.csv(saw.com.collect.alpine, "data/speciespersummit_miseqINFO.csv")

splist.phylo.info <- splist.raw.mi.tmp[!is.na(splist.raw.mi.tmp$Tip.Label.160511),] #424
dim(splist.phylo.info)

saw.com.mi <- spread(data = splist.phylo.info, key = Summit, value = Occurrence, fill = 0)
head(saw.com.mi)
dim(saw.com.mi)
saw.com.mi <- saw.com.mi[c(1, 10, 16:ncol(saw.com.mi))]
colnames(saw.com.mi) <- gsub(colnames(saw.com.mi), pattern = " ", replacement = "_")
write.csv(saw.com.mi, "data/speciespersummit_miseq.csv")


splist.phylo.alpine.info <- splist.raw.mi.tmp.alp[!is.na(splist.raw.mi.tmp.alp$Tip.Label.160511), c(1, 10, 11)] #424
dim(splist.phylo.alpine.info)  #307 speices total with MiSeq data
head(splist.phylo.alpine.info)

#### Talus
saw.com.alpine.mi <- spread(data = splist.phylo.alpine.info, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.alpine.mi)
colnames(saw.com.alpine.mi) <- gsub(colnames(saw.com.alpine.mi), pattern = " ", replacement = "_")

write.csv(saw.com.alpine.mi, "data/speciespersummit_miseqAlpine.csv")

#### Meadow
splist.phylo.meadow.info <- splist.raw.mi.tmp.mead[!is.na(splist.raw.mi.tmp.mead$Tip.Label.160511), c(1, 10, 11)] #424
dim(splist.phylo.meadow.info)  #117 speices total with MiSeq data
head(splist.phylo.meadow.info)

saw.com.meadow.mi <- spread(data = splist.phylo.meadow.info, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.meadow.mi)
colnames(saw.com.meadow.mi) <- gsub(colnames(saw.com.meadow.mi), pattern = " ", replacement = "_")

write.csv(saw.com.meadow.mi, "data/speciespersummit_miseqMeadow.csv")


############ Community Matrix MiSeq for Beta Diversity ############ 

head(splist.raw.mi.tmp)
splist.beta <- unique(splist.raw.mi.tmp[!is.na(splist.raw.mi.tmp$Annotated.Name), c("Annotated.Name", "Summit", "Occurrence")]) 
head(splist.beta)
saw.com.beta <- spread(data = splist.beta, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.beta) #143 unique species collected 
head(saw.com.beta)
colnames(saw.com.beta) <- gsub(colnames(saw.com.beta), pattern = " ", replacement = "_")
#write.csv(saw.com.beta, "data/speciespersummit_beta.csv")

## Talus
splist.beta.alp <- unique(splist.raw.mi.tmp[!is.na(splist.raw.mi.tmp$Annotated.Name), c("Annotated.Name", "Summit", "Meadow", "Occurrence")]) 
head(splist.beta.alp)
saw.com.beta.alp <- spread(data = splist.beta.alp %>% filter(Meadow==0), key = Summit, value = Occurrence, fill = 0)
dim(saw.com.beta.alp) #114 unique talus collected 
head(saw.com.beta.alp)
colnames(saw.com.beta.alp) <- gsub(colnames(saw.com.beta.alp), pattern = " ", replacement = "_")
#write.csv(saw.com.beta.alp[-2], "data/speciespersummit_beta_talus.csv")


############ Community Matrix PHLAWD ############ 
splist.phlawd.info <- splist.raw.mi.tmp[!splist.raw.mi.tmp$PHLAWD == "", c(5, 10, 11)] #313
dim(splist.phlawd.info) #313 species total from PHLAWD search
head(splist.phlawd.info)

splist <- splist.phlawd.info %>% group_by(Summit) %>% distinct(PHLAWD) 
saw.com.phlawd <- spread(data = splist, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.phlawd) #104

head(saw.com.phlawd)
colnames(saw.com.phlawd) <- gsub(colnames(saw.com.phlawd), pattern = " ", replacement = "_")
write.csv(saw.com.phlawd, "data/speciespersummit_phlawd.csv")



splist.phlawd.alpine.info <- splist.raw.mi.tmp.alp[!splist.raw.mi.tmp.alp$PHLAWD == "", c(4, 10, 11)] #313
dim(splist.phlawd.alpine.info) #313 species total from PHLAWD search
head(splist.phlawd.alpine.info)

splist <- splist.phlawd.alpine.info %>% group_by(Summit) %>% distinct(PHLAWD) 
saw.com.phlawd.alpine <- spread(data = splist, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.phlawd.alpine) #84

head(saw.com.phlawd.alpine)
colnames(saw.com.phlawd.alpine) <- gsub(colnames(saw.com.phlawd.alpine), pattern = " ", replacement = "_")

write.csv(saw.com.phlawd.alpine, "data/speciespersummit_phlawdAlpine.csv")

#### Meadow
splist.phlawd.meadow.info <- splist.raw.mi.tmp.mead[!splist.raw.mi.tmp.mead$PHLAWD == "", c(4, 10, 11)] #313
dim(splist.phlawd.meadow.info) #313 species total from PHLAWD search
head(splist.phlawd.meadow.info)

splist <- splist.phlawd.meadow.info %>% group_by(Summit) %>% distinct(PHLAWD) 
saw.com.phlawd.meadow <- spread(data = splist, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.phlawd.meadow) #57

head(saw.com.phlawd.meadow)
colnames(saw.com.phlawd.meadow) <- gsub(colnames(saw.com.phlawd.meadow), pattern = " ", replacement = "_")

write.csv(saw.com.phlawd.meadow, "data/speciespersummit_phlawdMeadow.csv")


############################################## 02_Clean ############################################## 
## for each alignment output for each gene region from PHLAWD, collapse intra-spacific taxa to species level, then keep longest sequence --> concatenate

source("R/ParsePHLAWD.R")
parsePHLAWD("output/01_PHLAWD/atpB/atpB.FINAL.aln.full") #8...8
parsePHLAWD("output/01_PHLAWD/ITS/ITS.FINAL.aln.full") #123...87
parsePHLAWD("output/01_PHLAWD/matK/matK.FINAL.aln.full") #71...62
parsePHLAWD("output/01_PHLAWD/rbcL/rbcL.FINAL.aln.full") #55..47
parsePHLAWD("output/01_PHLAWD/ndhF/ndhF.FINAL.aln.full") #19...19
parsePHLAWD("output/01_PHLAWD/trnTLF/trnTLF.FINAL.aln.full") #102...72


### Thanks https://uribeconvers.wordpress.com !! 
#### Make a list of alignment files in /02_Clean
# ls > names.txt

#Now make this file in TextWrangler easily by search and replace the following:
#  Search:(^.+)\r
#  Replace:phyutility -clean 0.5 -in \1 \2 -out \1\r
#  Search:$
#  Replace:.clean

# Now you have a file that you can make executable (chmod 777 names.txt) and run in the terminal (./names.txt)

# scp *.clean ../03_Concatenate/19genes/
# phyutility -concat -in *.clean -out sawtooth.phlawdALN.clean.concat.nex

### Convert from .nex to .fasta 
# NCLconverter sawtooth.phlawdALN.clean.concat.nex -efasta -osawtooth.phlawdALN.clean.concat
### Replace space in name
#sed -i "" 's/\ /_/g' sawtooth.phlawdALN.clean.concat.dna.fasta




############################################## 05_Trees ############################################## 
######## Use Congruifier (geiger) to ASSESS Trees 

genetree=read.nexus("output/05_Trees/concat/MiSeq/160511/sawtooth.160511.concat.nex") # saved as rooted in FigTree
tax=read.csv(file="output/06_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

tips=sapply(genetree$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""
atol=read.tree("output/06_Scaling/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
ftax=tax[match(atol$tip.label, rownames(tax)),]
ftax[,2]="Spermatophyta"
fatol=subset(atol, ftax, "family")
phy=genetree
swaptips=paste(1:length(tips),tips,sep="_")
phy$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
res
out=res$phy
congruif=out$node.label%in%res$calibrations$MRCA
out$node.label=NULL
out=nodelabel.phylo(out, tax, strict=FALSE)
out$node.label=ifelse(congruif, paste(out$node.label, "*", sep=""), out$node.label)
out$tip.label=genetree$tip.label[match(swaptips, res$phy$tip.label)] ### CHANGE tip labels back
out$tip.label=paste(genetree$tip.label[match(swaptips, res$phy$tip.label)], tax[,"family"], sep="=") ### ADD Family to tip labels

pdf("output/07_VisualizeTree/MiSeq/sawtooth.160511.concat.nex.congruifyTaxonomy.fan.pdf", width=10, height=10) 
tree <- ladderize(out, right=F)
plot.phylo(tree, type="fan", cex=0.05, label.offset = .05) 
nodelabels(out$node.label, frame="n", col="red", cex=0.2)
dev.off()

pdf("output/07_VisualizeTree/MiSeq/sawtooth.160511.concat.nex.congruifyTaxonomy.pdf") 
tree <- ladderize(out, right=F)
plot.phylo(edge.width = 0.25, tree, cex=0.05, label.offset = .05) 
nodelabels(out$node.label, frame="n", col="red", cex=0.2)
dev.off()

###help you evaluate where the tree is inconsistent with the taxonomy 
#see how the best node for lineage differs from the clade definition:
#out$FUN("Poaceae")
#missing from the clade in your tree OR unexpected within clade (but is found here in the tree)
out$FUN("Asteraceae")[[1]][[1]] #unexpected
"13_Unknown" "65_Unknown"

out$FUN("Orobanchaceae")

# * = Congruified Nodes
# Noded labels with ""  implies some inconsistency between the tips expected as defined in fleshed_genera.csv and the subtended tips in the tree at the best matching node


################### Label bootstrap support on unscaled ML tree
treeBS <- read.tree(file="output/05_Trees/concat/MiSeq/160511/Cipres_Data_sawtooth.160511.concat.MRE/RAxML_bipartitions.sawtooth.160511.concat.MRE") 
treeBS #424
## Find the node to root on
getMRCA(treeBS, tip=c("Pinus_albicaulis_14030", "Juniperus_communis_14330")) #647

treeBSroot <- root(treeBS, node=647, resolve.root=T)  
is.rooted(treeBSroot)
treeBSroot
#write.tree(treeBSroot, file="output/05_Trees/concat/MiSeq/160511/Cipres_Data_sawtooth.160511.concat.MRE/RAxML_bipartitions.sawtooth.160511.concat.MRE.unscaled.root.tre")

p1 <- character(length(treeBSroot$node.label)) # create a vector of colors that correspond to different support ranges
p1[treeBSroot$node.label >= 95] <- "turquoise3"
p1[treeBSroot$node.label < 95 & treeBSroot$node.label >= 75] <- "slategray"
p1[treeBSroot$node.label < 75] <- "red1"
p1[treeBSroot$node.label == 100] <- "green"
p1 
co = c("green", "turquoise3", "slategray", "red1")

pdf(file="output/07_VisualizeTree/MiSeq/Sawtooth.MiSeq.unscaled.BS.pdf") 
par(mar = rep(2, 4))
plot.phylo(ladderize(treeBSroot), type="fan", cex=0.2, label.offset=.01)
nodelabels(pch = 21, cex = 0.3, col=p1, bg = p1) 
legend("bottomleft", legend = expression(BS == 100, BS >= 95, 75 <= BS * " < 95", BS < 75), pch = 21, pt.bg = co, pt.cex = 1, cex=.5, bty = "n")
dev.off()


########################################################## 06_Scaling  ############################################################

##### WRITE treePL for TreeScaling: Re-run congruifier, changing scale = "NA", and writing out .treePL files

## Saved rooted in FigTree: rooted on all ferns & alies
genetree 

res1=congruify.phylo(fatol, phy, tax, tol=0, scale="NA") 

## WRITE OUT TREEPL FILE -- you'll need to do more than just run the exported file, but this gives you a start
nsites=8254 #SHOULD be number of nucleotides in alignment 
write.treePL(res1$target, res1$calibrations, nsites, base="output/06_Scaling/MiSeq/RAxML_bipartitions.sawtooth.160511.concat.MRE.treePL", opts=list(prime=TRUE))
## modify .infile to adjust location of .intree and outfile (local)
## Run on fortytwo: treePL EcrinSpPool.120515.1000.treePL.infile

#### After prime, changed end of block to estimated parameters
#PLACE THE LINES BELOW IN THE CONFIGURATION FILE
#opt = 5
#optad = 5
#optcvad = 3
#moredetailcvad


### Read in treePL output to change tip labels back
treePL <- read.tree("output/06_Scaling/MiSeq/RAxML_bipartitions.sawtooth.160511.concat.MRE.treePL.dated.tre")
treePL$tip.label=genetree$tip.label[match(swaptips, res1$target$tip.label)] ### CHANGE tip labels back

#Write .tre file
#write.tree(treePL, file="output/06_Scaling/MiSeq/RAxML_bipartitions.sawtooth.160511.concat.MRE.treePL.dated.rename.tre")


########### Match bootstrap node labels from ML tree to treePL scaled tree: Save as final Tree ###############
#Modified from: http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
source("R/treeFunctions.R")

treeBSroot #<- read.tree(file="output/05_Trees/concat/MiSeq/160511/Cipres_Data_sawtooth.160511.concat.MRE/RAxML_bipartitions.sawtooth.160511.concat.MRE.unscaled.root.tre) 
treePL #<- read.tree("output/06_Scaling/MiSeq/RAxML_bipartitions.sawtooth.160511.concat.MRE.treePL.dated.rename.tre")
treePL #1084

max(treePL$edge.length) # 265.5635

##### Match bootstrap node labels from ML tree to treePL scaled tree ####
## Takes a while for large trees...
plotTreePLBoot(treePL=treePL, bootTree=treeBSroot, file="output/06_Scaling/MiSeq/Sawtooth.160511.dated.bootstrap.tre") 
## Deleted "Root" label from nexus in TextWrangler 

################################################
########### Update Tip Names (Unknown sp., ...) and collapse to single species in Tree (need for beta-diversity)
## Note: branch lenghts are not euqal, because this is done on time-scaled tree
head(splist.raw.mi.tmp)
phy.beta.tmp <- sawMiseq$phy
phy.beta <- annotate.trim.tips(phy = phy.beta.tmp, spliton = "_", sepas = "_", taxonomy = splist.raw.mi.tmp)
phy.beta$tip.label
plot.phylo(phy.beta, cex=.2, type = "f")
write.tree(phy.beta, file = "data/sawtooth.miseq.beta.tre")


