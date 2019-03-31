
#### Plot community phylogeny with calibrated nodes, bootstrap support, 
#### pres/abs of trait data across tips, and the tip labels colored for introduced/native
#### Modified trait.plot function in diversitree package  (http://www.zoology.ubc.ca/~fitzjohn/diversitree.docs/trait.plot.html)

source("R/functions/diversitreeFunctions.R")
source("R/functions/treeFunctions.R")


plot.colorCom <- function (tree, dat, cols, lab = names(cols), col.names, class = NULL, 
                                 type = "f", w = w, legend = length(cols) > 1, cex.lab = 0.2, 
                                 font.lab = 3, cex.legend = 0.25, margin = 1/4, check = TRUE, 
                                 quiet = FALSE, cut.labs=NULL, leg.title = NULL, main = "", 
                                 leg.cex = 1, tip.labs = NULL, pch = 8, cex.shape = .25, ...) 
{
  if (!(type %in% c("f", "p"))) 
    stop("Only types 'f'an and 'p'hylogram are available")
  if (!is.null(class) && length(class) != length(tree$tip.label)) 
    stop("'class' must be a vector along tree$tip.label")
  n <- length(cols)
  if (n < 1) 
    stop("Need some colours")
  if (!all(tree$tip.label %in% rownames(dat))) 
    stop("All taxa must have entries in 'dat' (rownames)")
  if (n > 1) {
    if (!all(names(cols) %in% names(dat))) 
      stop("Not all colours have data")
    if (is.null(names(cols))) 
      stop("'cols' must be named")
    
  }

  dat <- dat[tree$tip.label, ,drop = FALSE]
  par(mar = rep(0, 4))
  t <- max(branching.times(tree))
  w <- w * t
  if (is.null(class)) {
    stopifnot(class(dat) == "data.frame", class(tree) == "phylo")
    len.tips <- length(tree$tip.label)
    len.taxa <- length(rownames(dat))
    if (len.tips != len.taxa | sum(tree$tip.label %in% rownames(dat)) != len.taxa) {
      stop("ERROR. Missing taxa in tree or data frame; # tips: ", 
           len.tips, "# taxa: ", len.taxa, "# tips in dat: ", 
           sum(tree$tip.label %in% rownames(dat)))
    }

    plt <- plot2.phylo(tree, type = type, show.tip.label = TRUE, 
                       label.offset = (n + 2) * w, cex = cex.lab, tip.color="black", ...)
  }
  if (type == "f") {
    xy <- plt$xy
    theta <- xy$theta[seq_along(tree$tip.label)]
    dt <- diff(sort(theta))[1]/2
    for (i in seq_along(cols)) {
        idx <- dat[[names(dat)[i]]]
        if (any(idx == 0, na.rm = TRUE)) 
          idx <- idx + 1
        filled.arcs(theta - dt, theta + dt, max(xy$x) + i * 
                      w, w, cols[[i]][idx]) 
      }
    }

  ######## Nodes Object ################
  #### Need to re-run Congruifier to get vec object of calibrated nodes
  vec <- get.cal.object(genetree=tree, taxonomy = "output/06_Scaling/Congruify/fleshed_genera.csv", refDates = "output/06_Scaling/Congruify/out_dates.tre")
  ##### Vector of BS supports 
  p2 <- character(length(tree$node.label)) # create a vector of colors that correspond to different support ranges 
  p2[] <- "#0000ff00" ## Transparent color: "#RRGGBBAA" and the AA portion is the opacity/trasparency.
  p2[tree$node.label >=  95] <- "black"
  p2[tree$node.label ==  100] <- "black"
  p2[tree$node.label < 95 & tree$node.label >= 75] <- "slategray4"
  p2[tree$node.label == ""] <- "black" # node = 353 = 100
  p2 
  
  nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
  nodelabels(pch = 21, cex = 0.3, col=p2, bg = p2) # Plot bootstrap support
  title(line = 0)
  if (is.null(cut.labs)) 
    cut.labs <- levs
  par(xpd=TRUE)
  #horiz=TRUE, pt.cex=1.5, bty = "n" inset = 0.025,
  
  invisible(plt)
}

#co <- c("black", "slategray4", "black", cols[[1]][2], cols[[2]][2], cols[[3]][2], cols[[4]][2], 
#cols[[5]][2], cols[[6]][2], cols[[7]][2], cols[[8]][2], cols[[9]][2], cols[[10]][2])


#legend(10,-1, legend = c(expression(BS >= 95, 75 <= BS * "< 95", Dated), 
#                         names(dat[2:ncol(dat)])), pch = c(21, 21, 1, rep(x=21, times= length(dat[2:ncol(dat)]))), pt.bg = co, 
#       cex=leg.cex, bty = "n", col=c("white", "white", "black", rep("white", times= length(dat[2:ncol(dat)])))) 
