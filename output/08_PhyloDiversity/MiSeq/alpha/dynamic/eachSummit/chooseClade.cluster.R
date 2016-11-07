##############################
## reduce communtiy phylogeny to clades of interest
tax=read.csv(file="data/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

###### Community Matrix Miseq for Beta phylo analyses with annotated species names and duplicate species removed (to the species level) ###### 
## All alpine (tallus + meadow)
sawtooth.com.miseq.beta <- read.csv(file="data/speciespersummit_beta.csv", row.names=2)
head(sawtooth.com.miseq.beta)
sawtooth.com.miseq.beta <- sawtooth.com.miseq.beta[2:ncol(sawtooth.com.miseq.beta)]
colnames(sawtooth.com.miseq.beta) <- gsub(colnames(sawtooth.com.miseq.beta) , pattern = "_", replacement = " ")
dim(sawtooth.com.miseq.beta) #143  unique species collected 


########## MiSeq Beta Phylogeny 
SawtoothMiseqBeta <- read.tree("data/sawtooth.miseq.beta.tre")
sawMiseqBeta <- comparative.comm(phy = SawtoothMiseqBeta, comm = as.matrix(t(sawtooth.com.miseq.beta))) 
sawMiseqBeta$comm

pruneCladeTaxonomyLookup <- function(tip.labels, tax, level, taxonomy){
  tips.ecrins=sapply(tip.labels, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips.ecrins, rownames(tax))
  ecrins_tax=tax[ll,]
  rownames(ecrins_tax)=names(tips.ecrins)
  ecrins_tax=as.matrix(ecrins_tax)
  ecrins_tax[is.na(ecrins_tax)]=""
  head(ecrins_tax)
  #length(which(ecrins_tax[,"Angiospermae"] == "Angiospermae")) # 1064 species are in Spermatophyta
  ecrins.clade <- names(which(ecrins_tax[,level] == taxonomy))
  #return(as.data.frame(ecrins_tax))
  return(ecrins.clade)
}

