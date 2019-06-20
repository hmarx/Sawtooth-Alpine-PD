#################################################################################################################
######## Preparing and cleaning datasets                      ###################################################
######## Hannah E. Marx, Jan 2019                             ###################################################
######## See workflow_20190331.rtf for more notes on pipeline ###################################################
#################################################################################################################


############# output/01_concatGeneSegments/
## Summarize the number of unique sequnces for species in MiSeq Dataset
parsePHLAWD("output/01_MiSeq/4_concensus/02_concateReads/atpB/atpB.miseq.aln.fasta") #246...107
parsePHLAWD("output/01_MiSeq/4_concensus/02_concateReads/ITS/ITS.R1_R2.aln.fasta") #187...100
parsePHLAWD("output/01_MiSeq/4_concensus/02_concateReads/matK/matK.miseq.aln.mod.fasta") #220...88
parsePHLAWD("output/01_MiSeq/4_concensus/02_concateReads/ndhF/ndhF.miseq.aln.fasta") #339...142
parsePHLAWD("output/01_MiSeq/4_concensus/02_concateReads/rbcL/rbcL.R1_R2.fasta") #85...62
parsePHLAWD("output/01_MiSeq/4_concensus/02_concateReads/trnTLF/trnTLF.miseq.aln.fasta") #136...81

# number of unique species in concatenated MiSeq Data
#  cat *unique.fasta > genesMiseq.fasta
parsePHLAWD(fasta.file = "output/01_MiSeq/4_concensus/02_concatGeneSegments/genesMiseq.fasta") #580...151

### Annotates .fasta alignment to updated nomenclature
source("R/functions/annotate.fasta.R") 
vouchertax <- read.csv("data/RAW/Sawtooth_CollectionVoucher20190420.csv") 
taxonomy.df <- cbind(as.character(vouchertax$Accepted.Name.GenBank), as.character(vouchertax$MiSeq.Label.190117))

rbcL.miseq <- annotate.fasta(fastaFile= "output/01_MiSeq/4_concensus/02_concatGeneSegments/rbcL.R1_R2.fasta", taxonomy = taxonomy.df) #62
atpB.miseq  <- annotate.fasta(fastaFile= "output/01_MiSeq/4_concensus/02_concatGeneSegments/atpB.miseq.aln.fasta", taxonomy = taxonomy.df) #107
matK.miseq  <- annotate.fasta(fastaFile= "output/01_MiSeq/4_concensus/02_concatGeneSegments/matK.miseq.aln.mod.fasta", taxonomy = taxonomy.df) #88
ndhF.miseq  <- annotate.fasta(fastaFile= "output/01_MiSeq/4_concensus/02_concatGeneSegments/ndhF.miseq.aln.fasta.rename.rem.fasta", taxonomy = taxonomy.df) #142
trnTLF.miseq  <- annotate.fasta(fastaFile= "output/01_MiSeq/4_concensus/02_concatGeneSegments/trnTLF.miseq.aln.fasta.rename.fasta", taxonomy = taxonomy.df) #81
ITS.miseq  <- annotate.fasta(fastaFile= "output/01_MiSeq/4_concensus/02_concatGeneSegments/ITS.R1_R2.aln.fasta.rename.fasta", taxonomy = taxonomy.df) #100

parsePHLAWD("output/01_MiSeq/4_concensus/02_concatGeneSegments/atpB.miseq.aln.fasta.rename.fasta") #246...107
parsePHLAWD("output/01_MiSeq/4_concensus/02_concatGeneSegments/matK.miseq.aln.mod.fasta.rename.fasta") #220...88
parsePHLAWD("output/01_MiSeq/4_concensus/02_concatGeneSegments/ITS.R1_R2.aln.fasta.rename.fasta.rename.fasta") #187...101
parsePHLAWD("output/01_MiSeq/4_concensus/02_concatGeneSegments/ndhF.miseq.aln.fasta.rename.rem.fasta.rename.fasta") #366...130
parsePHLAWD("output/01_MiSeq/4_concensus/02_concatGeneSegments/rbcL.R1_R2.fasta.rename.fasta") #85...62
parsePHLAWD("output/01_MiSeq/4_concensus/02_concatGeneSegments/trnTLF.miseq.aln.fasta.rename.fasta.rename.fasta") #136...82

# number of unique species in concatenated MiSeq Data
#  cat *unique.fasta > genesMiseq.renamed.fasta
parsePHLAWD(fasta.file = "output/01_MiSeq/4_concensus/02_concatGeneSegments/genesMiseq.renamed.fasta") #575...152 


############# output/02_GenBank/
## for each alignment output for each gene region from PHLAWD, collapse intra-spacific taxa to species level, then keep longest sequence --> concatenate
source("R/functions/ParsePHLAWD.R")
parsePHLAWD("output/02_GenBank/atpB/atpB.FINAL.aln.rn") #9...9
parsePHLAWD("output/02_GenBank/ITS/ITS.FINAL.aln.rn") #149...103
parsePHLAWD("output/02_GenBank/matK/matK.FINAL.aln.rn") #82...71 
parsePHLAWD("output/02_GenBank/rbcL/rbcL.FINAL.aln.rn") #55...47
parsePHLAWD("output/02_GenBank/ndhF/ndhF.FINAL.aln.rn") #24...23
parsePHLAWD("output/02_GenBank/trnTLF/trnTLF.FINAL.aln.rn") #118...84

# number of unique species in concatenated GenBank Data
#  cat *unique.fasta > genesGenBank.fasta
parsePHLAWD("output/02_GenBank/genesGenBank.fasta") #337...121


#################### Taxonomy updated to Voucher/GenBank (2019-04-01)
vouchertax <- read.csv("data/RAW/Sawtooth_CollectionVoucher20190420.csv") 
head(vouchertax)
#MiSeq_Label_190402 = updated
#MiSeq.Label.190117 = current

### Annotates .fasta alignment to updated nomenclature
source("R/functions/annotate.fasta.R") 

rbcL <- annotate.fasta.total(fastaFile= "output/06_geneTreesRemoved/rbcL/rbcL.total.fasta", taxonomy = vouchertax) #132
atpB <- annotate.fasta.total(fastaFile= "output/06_geneTreesRemoved/atpB/atpB.total.rem.fasta", taxonomy = vouchertax) #248
matK <- annotate.fasta.total(fastaFile= "output/06_geneTreesRemoved/matK/matK.total.rem.fasta", taxonomy = vouchertax) #271
ndhF <- annotate.fasta.total(fastaFile= "output/06_geneTreesRemoved/ndhF/ndhF.total.rem.fasta", taxonomy = vouchertax) #387
trnTLF <- annotate.fasta.total(fastaFile= "output/06_geneTreesRemoved/trnTLF/trnTLF.total.rem.fasta", taxonomy = vouchertax) #198
ITS <- annotate.fasta.total(fastaFile= "output/06_geneTreesRemoved/ITS/ITS.total.rem.fasta", taxonomy = vouchertax) #277

############# output/07_Concatenate/
# After cleaning, remove everythying in ID except Genus_species, and keep longest sequence:
source("R/functions/ParseTotal.R")
splist.raw <- read.csv("data/RAW/Sawtooth_CollectionVoucher20190420.csv") #first column already removed cf.
#splist.raw.sp <- gsub(pattern = "cf_", replacement = "", splist.raw[,1])

parsePHLAWDtotaldata(fasta.file = "output/07_Concatenate/atpB.total.rem.rename.aln.clean.fasta", splist.raw = as.character(splist.raw[,1])) #247...108
parsePHLAWDtotaldata(fasta.file = "output/07_Concatenate/ITS.total.rem.rename.aln.clean.fasta", splist.raw = splist.raw[,1]) #277...122
parsePHLAWDtotaldata(fasta.file = "output/07_Concatenate/matK.total.rem.rename.aln.clean.fasta", splist.raw = splist.raw[,1]) #271...111
parsePHLAWDtotaldata(fasta.file = "output/07_Concatenate/ndhF.total.rem.rename.aln.clean.fasta", splist.raw = splist.raw[,1]) #387...130
parsePHLAWDtotaldata(fasta.file = "output/07_Concatenate/rbcL.total.fasta.rename.aln.clean.fasta", splist.raw = splist.raw[,1]) #132...88
parsePHLAWDtotaldata(fasta.file = "output/07_Concatenate/trnTLF.total.rem.rename.aln.clean.fasta", splist.raw = splist.raw[,1]) #198...106


## save GenBank/MiSeq Identifiers for uniques sequences kept and used for each species  
source("R/functions/parseGBinfo.R")
atpbGB <- parseGBinfo(gene = "atpB", unique.GB.fasta = "output/07_Concatenate/atpB.total.unique.GB.fasta") #108
ITS <- parseGBinfo(gene = "ITS", unique.GB.fasta = "output/07_Concatenate/ITS.total.unique.GB.fasta") #122
tnrtlfGB <- parseGBinfo(gene = "trnTLF", unique.GB.fasta = "output/07_Concatenate/trnTLF.total.unique.GB.fasta") #106
rbclGB <- parseGBinfo(gene = "rbcL", unique.GB.fasta = "output/07_Concatenate/rbcL.total.unique.GB.fasta") #88
matkGB <- parseGBinfo(gene = "matK", unique.GB.fasta = "output/07_Concatenate/matK.total.unique.GB.fasta") #111
ndhfGB <- parseGBinfo(gene = "ndhF", unique.GB.fasta = "output/07_Concatenate/ndhF.total.unique.GB.fasta") #130

GB_accessions <- mergeGBinfo(ITS, tnrtlfGB, rbclGB, matkGB, atpbGB, ndhfGB)
dim(GB_accessions) #152 taxa, 6 gene regions
write.csv(GB_accessions, "output/07_Concatenate/GB_accessions.csv")


## Visualize taxonomy on tree:
genetree=read.nexus("output/08_Trees/sawtooth.totalData.clean.concat.tre") # saved as rooted in FigTree 
tax=read.csv(file="output/09_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)
tips=sapply(genetree$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""
atol=read.tree("output/09_Scaling/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
ftax=tax[match(atol$tip.label, rownames(tax)),]
ftax[,2]="Spermatophyta"
fatol=subset(atol, ftax, "family")
phy.total=genetree
swaptips=paste(1:length(tips),tips,sep="_")
phy.total$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res=congruify.phylo(fatol, phy.total, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
res
out=res$phy
congruif=out$node.label%in%res$calibrations$MRCA
out$node.label=NULL
out=nodelabel.phylo(out, tax, strict=FALSE)
out$node.label=ifelse(congruif, paste(out$node.label, "*", sep=""), out$node.label)
out$tip.label=genetree$tip.label[match(swaptips, res$phy$tip.label)] ### CHANGE tip labels back
out$tip.label=paste(genetree$tip.label[match(swaptips, res$phy$tip.label)], tax[,"family"], sep="=") ### ADD Family to tip labels

pdf("output/08_Trees/sawtooth.total.congruifyTaxonomy.fan.pdf", width=10, height=10) 
tree <- ladderize(out, right=F)
plot.phylo(tree, type="fan", cex=0.5, label.offset = .05) 
nodelabels(out$node.label, frame="n", col="red", cex=0.2)
dev.off()

pdf("output/08_Trees/sawtooth.total.congruifyTaxonomy.pdf") 
tree <- ladderize(out, right=F)
plot.phylo(edge.width = 0.25, tree, cex=0.05, label.offset = .05) 
nodelabels(out$node.label, frame="n", col="red", cex=0.2)
dev.off()

################### Label bootstrap support on unscaled ML tree
treeBStotal <- read.tree(file="output/08_Trees/CIPRES-sawtooth/RAxML_bipartitions.sawtooth.totalData.clean.concat") 
treeBStotal #152 taxa
## Find the node to root on
getMRCA(treeBStotal, tip=c("Selaginella_densa", "Selaginella_watsonii")) # 195

treeBStotalroot <- root(treeBStotal, node=195, resolve.root=T)  
is.rooted(treeBStotalroot)
treeBStotalroot
write.tree(treeBStotalroot, file="output/08_Trees/RAxML_bipartitions.sawtooth.total.unscaled.root.tre")

#### Drop ambiguous species:
#treeBStotalroot.tmp <- drop.tip(phy = treeBStotalroot, tip = "Erigeron")
#treeBStotalroot.tmp <- drop.tip(phy = treeBStotalroot.tmp, tip = "Phlox")
#treeBStotalroot.tmp <- drop.tip(phy = treeBStotalroot.tmp, tip = "Penstemon_sp")
#treeBStotalroot.tmp <- drop.tip(phy = treeBStotalroot.tmp, tip = "Unknown")
treeBStotalroot.tmp <- treeBStotalroot


### plot bootstrap support levels on phylo:
p1 <- character(length(treeBStotalroot.tmp$node.label)) # create a vector of colors that correspond to different support ranges
p1[treeBStotalroot.tmp$node.label >= 95] <- "turquoise3"
p1[treeBStotalroot.tmp$node.label < 95 & treeBStotalroot.tmp$node.label >= 75] <- "slategray"
p1[treeBStotalroot.tmp$node.label < 75] <- "red1"
p1[treeBStotalroot.tmp$node.label == 100] <- "green"
p1 
co = c("green", "turquoise3", "slategray", "red1")

pdf(file="output/08_Trees/Sawtooth.total.unscaled.BS.pdf") 
par(mar = rep(2, 4))
plot.phylo(ladderize(treeBStotalroot.tmp), type="fan", cex=0.2, label.offset=.01)
nodelabels(pch = 21, cex = 0.3, col=p1, bg = p1) 
legend("bottomleft", legend = expression(BS == 100, BS >= 95, 75 <= BS * " < 95", BS < 75), pch = 21, pt.bg = co, pt.cex = 1, cex=.5, bty = "n")
dev.off()


################## output/09_Scaling/
##### WRITE treePL for TreeScaling: Re-run congruifier, changing scale = "NA", and writing out .treePL files
res1.total=congruify.phylo(fatol, phy.total, tax, tol=0, scale="NA") 

## WRITE OUT TREEPL FILE -- you'll need to do more than just run the exported file, but this gives you a start
nsites=3193 #SHOULD be number of nucleotides in alignment 
write.treePL(res1.total$target, res1.total$calibrations, nsites, base="output/09_Scaling/RAxML_bipartitions.sawtooth.concat.treePL", opts=list(prime=TRUE))
## modify .infile to adjust location of .intree and outfile (local)
## $ treePL RAxML_bipartitions.sawtooth.concat.treePL.infile

#### After prime, changed end of block to estimated parameters
#PLACE THE LINES BELOW IN THE CONFIGURATION FILE
#opt = 2
#optad = 2
#optcvad = 2

### Read in treePL output to change tip labels back
treePL.total <- read.tree("output/09_Scaling/RAxML_bipartitions.sawtooth.concat.treePL.dated.tre")
treePL.total$tip.label=genetree$tip.label[match(swaptips, res1.total$target$tip.label)] ### CHANGE tip labels back

#Write .tre file
write.tree(treePL.total, file="output/09_Scaling/RAxML_bipartitions.sawtooth.concat.treePL.dated.rename.tre")

########### Match bootstrap node labels from ML tree to treePL scaled tree: Save as final Tree!! ###############
#Modified from: http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
source("R/functions/treeFunctions.R")

treeBSroot.total <- read.tree(file="output/08_Trees/RAxML_bipartitions.sawtooth.total.unscaled.root.tre")
treePL.total #<- read.tree("output/09_Scaling/RAxML_bipartitions.sawtooth.concat.treePL.dated.rename.tre")
treePL.total # 152
max(treePL.total$edge.length) # 550.9611

##### Match bootstrap node labels from ML tree to treePL scaled tree ####
## Takes a while for large trees...
plotTreePLBoot(treePL=treePL.total, bootTree=treeBSroot.total, file="output/09_Scaling/Sawtooth.total.dated.bootstrap.tre") 
## Deleted "Root" label from nexus in TextWrangler 

sawtooth.tree <- read.tree("output/09_Scaling/Sawtooth.total.dated.bootstrap.tre")
#### Drop ambiguous species:
sawtooth.tree.tmp <- drop.tip(phy = sawtooth.tree, tip = "Erigeron")
sawtooth.tree.tmp <- drop.tip(phy = sawtooth.tree.tmp, tip = "Phlox")
sawtooth.tree.tmp <- drop.tip(phy = sawtooth.tree.tmp, tip = "Draba_crassifolia") ## This was recenlty annotated as ambiguous

#sawtooth.tree.tmp #149 taxa

# Write Final .tre file!! 
write.tree(sawtooth.tree.tmp, file="output/09_Scaling/sawtooth.totalData.tre")




################################################## Community Matrix ################################################## 
splist.raw.mi <- read.csv("data/RAW/Sawtooth_CollectionVoucher20190420.csv", as.is=T) ### Updated with Genbank Taxonomy 02 Feb 2019

head(splist.raw.mi)
dim(splist.raw.mi) # 476 species collected

length(unique(na.omit(splist.raw.mi.tmp$MiSeq.Accession.Number))) # 419 MiSeq accessions

splist.raw.mi.tmp <- filter(splist.raw.mi, Occurrence == 0)
dim(splist.raw.mi.tmp) # 476

splist.raw.mi.tmp.alp <- splist.raw.mi.tmp[splist.raw.mi.tmp$Meadow != 1,] #true alpine species 
dim(splist.raw.mi.tmp.alp) # 359 talus species
splist.raw.mi.tmp.mead <- splist.raw.mi.tmp[splist.raw.mi.tmp$Meadow == 1,] #meadow species
dim(splist.raw.mi.tmp.mead) # 117 meadow species 


########################## Summarize Collections, MiSeq and PHLAWD
head(splist.raw.mi.tmp)
sum1 <- splist.raw.mi.tmp %>% group_by(Summit) %>% dplyr::summarise(n_collected = n_distinct(Collection), 
                                                                    n_MiSeq_individuals = n_distinct(MiSeq.Accession.Number), 
                                                                    n_MiSeq_species = n_distinct(Accepted.Name.GenBank.MiSeq.Data),
                                                                    n_PHLAWD = n_distinct(PHLAWD))

sum2 <- filter(splist.raw.mi.tmp, Meadow == 0) %>%  group_by(Summit) %>% dplyr::summarise(n_talus_collected = n_distinct(Collection), 
                                                                                          n_talus_individuals = n_distinct(MiSeq.Accession.Number), 
                                                                                          n_talus_species = n_distinct(Accepted.Name.GenBank),
                                                                                          n_talus_PHLAWD = n_distinct(PHLAWD))

sum3 <- filter(splist.raw.mi.tmp, Meadow == 1) %>%  group_by(Summit) %>% dplyr::summarise(n_meadow_collected = n_distinct(Collection), 
                                                                                          n_meadow_individuals = n_distinct(MiSeq.Accession.Number), 
                                                                                          n_meadow_species = n_distinct(Accepted.Name.GenBank),
                                                                                          n_meadow_PHLAWD = n_distinct(PHLAWD))

sum4 <- splist.raw.mi.tmp %>% dplyr::summarise(n_collected = n_distinct(Collection), n_MiSeq_individuals = n_distinct(na.omit(MiSeq.Accession.Number)), 
                                               n_MiSeq_species = n_distinct(Accepted.Name.GenBank.MiSeq.Data),
                                               n_PHLAWD = n_distinct(PHLAWD))

sum5 <- filter(splist.raw.mi.tmp, Meadow == 0) %>%  dplyr::summarise(n_talus_collected = n_distinct(Collection), 
                                                                     n_talus_individuals = n_distinct(na.omit(MiSeq.Accession.Number)), 
                                                                     n_talus_species = n_distinct(Accepted.Name.GenBank),
                                                                     n_talus_PHLAWD = n_distinct(PHLAWD))

sum6 <- filter(splist.raw.mi.tmp, Meadow == 1) %>%  plyr::summarise(n_meadow_collected = n_distinct(Collection), 
                                                                    n_meadow_individuals = n_distinct(na.omit(MiSeq.Accession.Number)), 
                                                                    n_meadow_species = n_distinct(Accepted.Name.GenBank),
                                                                    n_meadow_PHLAWD = n_distinct(PHLAWD))

sumtot <- merge(cbind(sum1, sum2[-1]), sum3, by =1, all.x=T)
sumtot <- rbind(sumtot, cbind(Summit = "Total", sum4, sum5,sum6))
write.csv(sumtot, file="figs/tables/summary.collection.csv")



############ Community Total Data ############ 


#### All Species collected
head(splist.raw.mi.tmp)
splist.raw.mi.tmp$Accepted.Name.GenBank <- gsub(pattern = "cf_", replacement = "", x = splist.raw.mi.tmp$Accepted.Name.GenBank)
splist.col <- splist.raw.mi.tmp %>% group_by(Summit, Occurrence, Meadow) %>% distinct(Accepted.Name.GenBank) 
saw.com.collect.all <- spread(data = splist.col, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.collect.all) # 199 total species collected
saw.com.collect.all <- saw.com.collect.all %>% group_by(Accepted.Name.GenBank) %>% dplyr::summarise_all(funs(max))
write.csv(saw.com.collect.all, "data/speciespersummit_collect.csv")
write.csv(merge(GB_accessions, saw.com.collect.all, by.x=0, by.y=1, all=T), file = "figs/tables/accession.collections.csv")


#### Talus
saw.com.collect.alpine.tmp <- splist.col[splist.col$Meadow == 0,]
saw.com.collect.alpine <- spread(data = saw.com.collect.alpine.tmp, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.collect.alpine) # 123 alpine species collected 
head(saw.com.collect.alpine)
write.csv(saw.com.collect.alpine, "data/speciespersummit_collectAlpine.csv")

#### Meadow
splist.collect.meadow.info <- splist.col[splist.col$Meadow ==1,]
head(splist.collect.meadow.info)
saw.com.collect.meadow <- spread(data = splist.collect.meadow.info, key = Summit, value = Occurrence, fill = 0)
dim(saw.com.collect.meadow) # 76 meadow species collected 
head(saw.com.collect.meadow)
write.csv(saw.com.collect.meadow, "data/speciespersummit_collectMeadow.csv")


