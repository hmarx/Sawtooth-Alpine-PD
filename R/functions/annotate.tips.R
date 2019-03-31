## Annotates tip names 

#phy=phy.tmp
#spliton = "_"
#sepas = " "
#taxonomy = df <- c("updated", "tips")
### Get Just Genus_species...(remove infraspecific identifiers) 
annotate.tips <- function(phy, spliton, sepas, taxonomy){
  split <- strsplit(as.character(phy$tip.label), split=spliton, fixed=TRUE)
  genus.name <- sapply(split, "[", 1L) #get just the genus_species_var...
  species.name <- sapply(split, "[", 2L) #get just the species_var...
  combinedname <- paste(genus.name, species.name, sep=sepas) #get just genus_species
  tab <- as.data.frame(cbind(tips = phy$tip.label, combinedname))
  ## Make a taxonomy lookup list 
  tmp.merge <- merge(tab, taxonomy[c(1:2)], by.x= 1, by.y=2, all.x=T, all.y=F)
  
  for (i in 1:length(phy$tip.label)){
    print(i)
    if (phy$tip.label[i] %in% tmp.merge$tips){
      x <- tmp.merge[which(tmp.merge$tips==phy$tip.label[i]), "tips"]
      xx <- tmp.merge[which(tmp.merge$tips==phy$tip.label[i]), "combinedname"]
      y <- tmp.merge[which(tmp.merge$tips==phy$tip.label[i]), 3]
      phy$tip.label[i] <- y
      if (xx != y){
        print(paste(i, "input =", x, ", accepted =", y))
        
      }
      
    }
  }
  
  # get just genus species 
  #phy$tip.label<-sapply(strsplit(phy$tip.label,"_"),function(x) paste(x[1], x[2], sep="_"))
  
  return(phy) 
  
}

### The same as above but also removes duplicate species
annotate.trim.tips <- function(phy, spliton, sepas, taxonomy){
  split <- strsplit(as.character(phy$tip.label), split=spliton, fixed=TRUE)
  genus.name <- sapply(split, "[", 1L) #get just the genus_species_var...
  species.name <- sapply(split, "[", 2L) #get just the species_var...
  combinedname <- paste(genus.name, species.name, sep=sepas) #get just genus_species
  tab <- as.data.frame(cbind(tips = phy$tip.label, combinedname))
  ## Make a taxonomy lookup list 
  tmp.merge <- merge(tab, taxonomy[c(1:2)], by.x= 1, by.y=2, all.x=T, all.y=F)
  ## here are our unique species in the alpine communies
  tmp.merge <- tmp.merge[!duplicated(tmp.merge$Annotated_Name),]
  head(tmp.merge)
  ## now drop all but one of each
  phy.prune <- drop.tip(phy = phy, tip = setdiff(phy$tip.label, tmp.merge[[1]]))

  for (i in 1:length(phy.prune$tip.label)){
    print(i)
    if (phy.prune$tip.label[i] %in% tmp.merge$tips){
      x <- tmp.merge[which(tmp.merge$tips==phy.prune$tip.label[i]), "tips"]
      xx <- tmp.merge[which(tmp.merge$tips==phy.prune$tip.label[i]), "combinedname"]
      y <- tmp.merge[which(tmp.merge$tips==phy.prune$tip.label[i]), 3]
      phy.prune$tip.label[i] <- y
      if (xx != y){
        print(paste(i, "input =", x, ", accepted =", y))
        
      }
      
    }
  }

  # get just genus species 
  #phy.prune$tip.label<-sapply(strsplit(phy.prune$tip.label,"_"),function(x) paste(x[1], x[2], sep="_"))
  
  return(phy.prune) 
  
}



