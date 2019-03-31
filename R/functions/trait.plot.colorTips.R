
#### Plot community phylogeny with calibrated nodes, bootstrap support, 
#### pres/abs of trait data across tips, and the tip labels colored for introduced/native
#### Modified trait.plot function in diversitree package  (http://www.zoology.ubc.ca/~fitzjohn/diversitree.docs/trait.plot.html)

source("R/functions/diversitreeFunctions.R")
source("R/functions/treeFunctions.R")


trait.plot.colorTip <- function (tree, dat, cols, shape1, shape2, pch1, pch2, lab = names(cols), datTr, trait, taxa.names, str = NULL, class = NULL, 
          type = "f", w = 1/50, legend = length(cols) > 1, cex.lab = 0.5, 
          font.lab = 3, cex.legend = 0.75, margin = 1/4, check = TRUE, 
          quiet = FALSE, cut.labs=NULL, leg.title = NULL, main = "", 
          leg.cex = 1, tip.labs = NULL,cex.shape, ...) 
{
  if (!(type %in% c("f", "p"))) 
    stop("Only types 'f'an and 'p'hylogram are available")
  if (!is.null(class) && length(class) != length(tree$tip.label)) 
    stop("'class' must be a vector along tree$tip.label")
  n <- length(cols)
  if (n < 1) 
    stop("Need some colours")
  if (!is.data.frame(dat)) {
    if (is.vector(dat) && n == 1) {
      nm <- names(dat)
      dat <- matrix(dat)
      rownames(dat) <- nm
    }
    else {
      stop("dat must be a matrix")
    }
  }
  if (!all(tree$tip.label %in% rownames(dat))) 
    stop("All taxa must have entries in 'dat' (rownames)")
  if (n > 1) {
    if (!all(names(cols) %in% names(dat))) 
      stop("Not all colours have data")
    if (is.null(names(cols))) 
      stop("'cols' must be named")
    dat <- dat[names(cols)]
  }
  if (is.null(str)) {
    str <- lapply(dat, function(x) as.character(sort(unique(x))))
  }
  dat <- dat[tree$tip.label, , drop = FALSE]
  par(mar = rep(0, 4))
  t <- max(branching.times(tree))
  w <- w * t
  if (is.null(class)) {
    stopifnot(trait %in% names(datTr), taxa.names %in% names(datTr), 
              class(datTr) == "data.frame", class(tree) == "phylo")
    len.tips <- length(tree$tip.label)
    len.taxa <- length(datTr[, taxa.names])
    if (len.tips != len.taxa | sum(tree$tip.label %in% datTr[,taxa.names]) != len.taxa) {
      stop("ERROR. Missing taxa in tree or data frame; # tips: ", 
           len.tips, "# taxa: ", len.taxa, "# tips in datTr: ", 
           sum(tree$tip.label %in% datTr[, taxa.names]))
    }
    order <- match(tree$tip.label, datTr[, taxa.names])
    ordered.trait <- datTr[trait][order, ]
    if (is.factor(ordered.trait)) {
      levs <- levels(ordered.trait)
      tip.color <- rep("black", times = len.taxa)
    }
    else {
      tip.color = "black"
    }
    #if (!is.null(tip.labs)) {
    #  tree$tip.label <- datTr[tip.labs][order, ]
    #}
    
    plt <- plot2.phylo(tree, type = type, show.tip.label = TRUE, 
                       label.offset = (n + 3) * w, cex = cex.lab, tip.color=tip.color, ...)
  }
  else {
    plt <- plot2.phylo(tree, type = type, show.tip.label = FALSE, 
                       label.offset = t * margin, ...)
    group.label.tip(plt, class, "black", "black", lwd = 1.5, 
                    offset.bar = w * (n + 2), offset.lab = w * (n + 4), 
                    cex = cex.lab, font = font.lab, check = check, quiet = quiet)
  }
  if (type == "f") {
    xy <- plt$xy
    theta <- xy$theta[seq_along(tree$tip.label)]
    dt <- diff(sort(theta))[1]/2
    for (i in seq_along(cols)) {
      if (names(cols[i]) %in% shape1){
        tip = which(t(dat)[names(dat)[i],] > 0)
        col = cols[[i]][[2]]
        #lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        XX <- plt$xx[tip]
        YY <- plt$yy[tip]
        BOTHlabels(XX = XX * 1.24, YY = YY * 1.24, adj = c(0,0), frame = "rect",  
                   thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
                   bg = "black", horiz = FALSE, width = FALSE, height = NULL, pch = pch1, cex = cex.shape)
      } else if (names(cols[i]) %in% shape2){
        tip = which(t(dat)[names(dat)[i],] > 0)
        col = cols[[i]][[2]]
        #lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        XX <- plt$xx[tip]
        YY <- plt$yy[tip]
        BOTHlabels(XX = XX * 1.27, YY = YY * 1.27, adj = c(0,0), frame = "rect",  
                   thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
                   bg = "black", horiz = FALSE, width = FALSE, height = NULL, pch = pch2, cex = cex.shape)
      } else{
        idx <- dat[[names(dat)[i]]]
        if (any(idx == 0, na.rm = TRUE)) 
          idx <- idx + 1
        filled.arcs(theta - dt, theta + dt, max(xy$x) + i * 
                      w, w, cols[[i]][idx])
        
      }

    }
  }
  else {
    xy <- plt$xy[seq_along(tree$tip.label), ]
    dy <- 0.5
    for (i in seq_along(cols)) {
      idx <- dat[[names(dat)[i]]]
      if (any(idx == 0, na.rm = TRUE)) 
        idx <- idx + 1
      xleft <- xy[1, 1] + w * i
      xright <- xleft + w
      ybottom <- xy[, 2] - dy
      ytop <- ybottom + dy * 2
      rect(xleft, ybottom, xright, ytop, col = cols[[i]][idx], 
           border = NA)
    }
  }
  ######## Nodes Object ################
  #### Need to re-run Congruifier to get vec object of calibrated nodes
  vec <- get.cal.object(genetree=tree, taxonomy = taxonomyDir, refDates = refDates)
  ##### Vector of BS supports 
  p2 <- character(length(tree$node.label)) # create a vector of colors that correspond to different support ranges 
  p2[] <- "#0000ff00" ## Transparent color: "#RRGGBBAA" and the AA portion is the opacity/trasparency.
  p2[tree$node.label >=  95] <- "black"
  p2[tree$node.label ==  100] <- "black"
  p2[tree$node.label < 95 & tree$node.label >= 75] <- "slategray4"
  p2[tree$node.label == ""] <- "black" # node = 353 = 100
  p2 
  
  nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
  nodelabels(pch = 21, cex = 0.5, col=p2, bg = p2) # Plot bootstrap support
  title(line = 0)
  if (is.null(cut.labs)) 
    cut.labs <- levs

  invisible(plt)
}

#co <- c("black", "slategray4", "black", cols[[1]][1], cols[[1]][2], cols[[2]][2], cols[[3]][2], cols[[4]][2], cols[[5]][2], cols[[6]][2])
#legend("bottomright", cut.labs, legend = c(expression(BS >= 95, 75 <= BS * " < 95", Dated), 
#                                          "Invasive", "Native","Seed Mass", "Max Height",  "SLA", "Leaf Size", "Leaf N"), pch = c(21, 21, 1, 21,21,21,21,21), pt.bg = co, 
#       inset = 0.0025, cex=leg.cex, bty = "n", col=c("white", "white", "black", "white","white","white","white","white")) #pt.cex=1.5, bty = "n"
