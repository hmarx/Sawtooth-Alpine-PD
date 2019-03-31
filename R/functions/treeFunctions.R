

##The getAllSubTrees function below is a necessary subfunction that atomizes a tree into each individual subclade and was provided compliments of Luke Harmon.
getAllSubtrees<-function(phy, minSize=2) {
  res<-list()
  count=1
  ntip<-length(phy$tip.label)
  for(i in 1:phy$Nnode) {
    l<-tips(phy, ntip+i)
    bt<-match(phy$tip.label, l)
    if(sum(is.na(bt))==0) {
      st<-phy} 
    else st<-drop.tip(phy, phy$tip.label[is.na(bt)])
    if(length(st$tip.label)>=minSize) {
      res[[count]]<-st
      count<-count+1
    }
  }
  res
}

plotTreePLBoot <- function(treePL,bootTree, file) {
  getAllSubtrees(treePL)->treePLSub
  getAllSubtrees(bootTree)->bootSub
  #bootList<-matrix("<50",Nnode(treePL),1)
  bootList<-matrix("", Nnode(treePL),1)
  
  #The commands below compare all the subclades in the Bayes tree to all the subclades in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
  for(i in 1:Nnode(treePL)) {
    for(j in 1:Nnode(bootTree)) {
      match(treePLSub[[i]]$tip.label[order(treePLSub[[i]]$tip.label)], bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)])->shared
      match(bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)], treePLSub[[i]]$tip.label[order(treePLSub[[i]]$tip.label)])->shared2
      if(sum(is.na(c(shared,shared2)))==0) {
        bootTree$node.label[j]->bootList[i]
      }}}
  treePLBS <- treePL
  treePLBS$node.label <- bootList
  plot(ladderize(treePLBS, right=F), cex=.10, lwd=0.1) #Plots your Bayesian consensus tree
  nodelabels(treePLBS$node.label, adj=c(1.2, -0.3), frame="n", cex=.2, font=2) #Adds bootstrap values.
  tree <- as.phylo(treePLBS)
  write.tree(tree, file=file) #SAVE as tree
  
}


#### Vector of congruified nodes
mrcaID=function(phy, cal){
  cal=as.matrix(cal)
  res=sapply(1:nrow(cal), function(idx){ ## loop over rows
    tips=cal[idx, c("taxonA", "taxonB")] ## fetch spanning taxa
    return(geiger:::.mrca(tips, phy)) ## MRCA of spanning taxa (node ID)
  })
  N=Ntip(phy)
  n=Nnode(phy)
  nn=integer(N+n) ## create empty vector of same length as branches in tree
  nn[res]=1 ## identify nodes that appear within calibrations
  nn=nn[-c(1:N)] ## exclude tip branches
  return(nn) ## return vector (ordered from first to last internal node in the tree)
}




#### Without Traits
#pdf("SJtreePLbootsStatusDropNEW.pdf")
#color.plot.phylo3NoTrait(treePLboots, tra, "status", "Row.names", col.names = c("magenta1", "green4"), label.offset=.02, phy.cex=.3, cut.labs =c("Introduced", "Native"), leg.cex=.8) ##"turquoise4", "tomato"
#dev.off()

color.plot.phylo3NoTrait <- function (phylo, df, trait, taxa.names, label.offset, phy.cex=.3, num.breaks = ifelse(is.factor(df[, trait]), length(levels(df[, trait])), 12), col.names = rainbow(ifelse(length(num.breaks) > 1, length(num.breaks) - 1, num.breaks)), cut.labs = NULL, leg.title = NULL, main = trait, leg.cex = 1, tip.labs = NULL, ...) 
{
  init.par <- par(mar = c(0, 0, 1, 0))
  stopifnot(trait %in% names(df), taxa.names %in% names(df), 
            class(df) == "data.frame", class(phylo) == "phylo")
  len.tips <- length(phylo$tip.label)
  len.taxa <- length(df[, taxa.names])
  if (len.tips != len.taxa | sum(phylo$tip.label %in% df[,taxa.names]) != len.taxa) {
    stop("ERROR. Missing taxa in tree or data frame; # tips: ", 
         len.tips, "# taxa: ", len.taxa, "# tips in df: ", 
         sum(phylo$tip.label %in% df[, taxa.names]))
  }
  order <- match(phylo$tip.label, df[, taxa.names])
  ordered.trait <- df[trait][order, ]
  if (is.factor(ordered.trait)) {
    levs <- levels(ordered.trait)
    tip.color <- rep("black", times = len.taxa)
    tip.color <- col.names[match(ordered.trait, levs)]
  }
  else {
    tip.color = as.character(cut(ordered.trait, breaks = num.breaks, 
                                 labels = col.names))
    levs <- levels(cut(ordered.trait, breaks = num.breaks))
  }
  if (!is.null(tip.labs)) {
    phylo$tip.label <- df[tip.labs][order, ]
  }
  
  
  plot.phylo(ladderize(phylo), type="fan", cex=phy.cex, label.offset=label.offset, tip.color=tip.color, main=main,...) #label.offset=.02
  nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
  nodelabels(pch = 21, cex = 0.3, col=p2, bg = p2) # Plot bootstrap support
  title(line = 0)
  if (is.null(cut.labs)) 
    cut.labs <- levs
  legend("topleft", cut.labs, fill = col.names, inset = 0.025, 
         title = leg.title, cex = leg.cex, bty = "n")
  co <- c("black", "slategray4", "black")
  legend("bottomright", cut.labs, legend = c(expression(BS >= 95, 75 <= BS * " < 95", Dated)), pch = c(21, 21, 1, 21,21,21,21,21), pt.bg = co, 
         inset = 0.0025, cex=leg.cex, bty = "n", col=c("white", "white", "black", "white")) #pt.cex=1.5, bty = "n"
  #co2 <- c("orange", "purple", "gold1", "green4", "red2")
  #legend("bottomleft", cut.labs, legend = c("Leaf N", "SLA", "Max Height", "Leaf Size", "Seed Mass"), pch = c(21, 21, 1), pt.bg = co2, 
  #inset = 0.025, cex=leg.cex) #pt.cex=1.5, bty = "n"
  on.exit(par(init.par))
}



get.cal.object <- function (genetree, taxonomy, refDates){
  genetree=genetree
  tax=read.csv(file=taxonomy, as.is=TRUE, row=1) 
  head(tax)
  tips=sapply(genetree$tip.label, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips, rownames(tax))
  SJ_tax=tax[ll,]
  head(SJ_tax)
  rownames(SJ_tax)=names(tips)
  SJ_tax=as.matrix(SJ_tax)
  SJ_tax[is.na(SJ_tax)]=""
  atol=read.tree(refDates) #dataed reference tree, Soltis et al. 2011
  ftax=tax[match(atol$tip.label, rownames(tax)),]
  head(ftax)
  ftax[,2]="Spermatophyta"
  fatol=subset(atol, ftax, "family")
  phy=genetree
  swaptips=paste(1:length(tips),tips,sep="_")
  phy$tip.label=swaptips
  tax=SJ_tax
  rownames(tax)=swaptips
  res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
  res
  cal=res$calibrations
  #### Vector of congruified nodes
  vec=mrcaID(phy, cal) ## get vector for calibrated nodes
  sum(vec)==nrow(cal) ## check on whether the function is working appropriately
  # plot.phylo(phy, type="fan", show.tip=FALSE, edge.width=0.1)
  ## plot box at node only if calibrated
  # nodelabels(text=NULL, cex=ifelse(vec==0, NA, 2), frame="n", bg="lightskyblue", col="lightgray", pch=22)
  return(vec)

  
}

