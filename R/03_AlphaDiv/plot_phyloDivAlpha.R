#####################################################################################################################
############# Plots of phylogenetic diversity within alpine summits (alpha) ################################################## 
############# Hannah E. Marx, 20 Jan 2019 ########################################################################### 
#####################################################################################################################

############## Tile plots of alpha SES divided by clade

master.ses.alpha <- read.csv(file="output/10_PhyloDiversity/alpha/sawtooth.alpha.RD.SES.csv", row.names=1)
head(master.ses.alpha)
str(master.ses.alpha)
master.ses.alpha$sig
master.ses.alpha$metric <- factor(master.ses.alpha$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))
master.ses.alpha$data <- factor(master.ses.alpha$data, labels = c("total"="Total Data"))
master.ses.alpha$clade <- factor(master.ses.alpha$clade, levels = c("Ericales", "Brassicales", "Lamiales", "Caryophyllales", "Poales", "Asterales", "Tracheophyta"))
master.ses.alpha$pool <- factor(master.ses.alpha$pool, levels = c("All Alpine", "Talus", "Meadow"))
master.ses.alpha$summits <- factor(master.ses.alpha$summits, levels = c( "Horstmann Peak", "Braxon Peak","Thompson Peak", "Snowyside Peak", "Mount Cramer",  "D.O. Lee Peak","Salzburger Spitzl",  "Castle Peak",  "Hyndman Peak"))
master.ses.alpha$data <- factor(master.ses.alpha$data, levels = c("Total Data"))

sort(master.ses.alpha$obs.z) #-3.194197455, 2.231214848

### MNTD
plot.pools <- ggplot(filter(master.ses.alpha, metric == "MNTD"), aes(y=clade, x=summits, fill=as.numeric(as.character(obs.z))))
plot.pools <- plot.pools + geom_tile(colour="white", alpha=.75)
plot.pools <- plot.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-4, 4))
plot.pools <- plot.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.pools <- plot.pools + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.pools <- plot.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.pools <- plot.pools + facet_grid(pool ~ data)#,space="free",scales="free", as.table = F)   
plot.pools <- plot.pools + scale_x_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + scale_y_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = 10 * 0.8, angle = -45, hjust = 0, colour = "black"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="figs/alpha/alpha_SES_MNTD_tile.pdf")
plot.pools
dev.off()

### MPD
plot.pools <- ggplot(filter(master.ses.alpha, metric == "MPD"), aes(y=clade, x=summits, fill=obs.z))
plot.pools <- plot.pools + geom_tile(colour="white", alpha=.75)
plot.pools <- plot.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-4, 4))
plot.pools <- plot.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.pools <- plot.pools + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.pools <- plot.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.pools <- plot.pools + facet_grid(pool ~ data)#,space="free",scales="free", as.table = F)   
plot.pools <- plot.pools + scale_x_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + scale_y_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = 10 * 0.8, angle = -45, hjust = 0, colour = "black"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="figs/alpha/alpha_SES_MPD_tile.pdf")
plot.pools
dev.off()


############## Species in All Apline
plot.pools <- ggplot(filter(master.ses.alpha, data == "Total Data" & pool == "All Alpine"), aes(y=clade, x=summits, fill=obs.z))
plot.pools <- plot.pools + geom_tile(colour="white", alpha=.75)
plot.pools <- plot.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-4, 4))
plot.pools <- plot.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.pools <- plot.pools + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.pools <- plot.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.pools <- plot.pools + facet_grid(metric ~ .)#,space="free",scales="free", as.table = F)   
plot.pools <- plot.pools + scale_x_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + scale_y_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = 10 * 0.8, angle = -45, hjust = 0, colour = "black"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="figs/alpha/alpha_SES_TotalData_tile.pdf")
plot.pools
dev.off()


