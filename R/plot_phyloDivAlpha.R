############## Tile plots of alpha SES divided by clade

master.ses.alpha <- read.csv(file="output/08_PhyloDiversity/MiSeq/alpha/static/alpine.phylogeny.pool.SES.csv", row.names=1)
head(master.ses.alpha)
str(master.ses.alpha)
master.ses.alpha$sig
master.ses.alpha$clade <- factor(master.ses.alpha$clade, levels = c("Ericales", "Brassicales", "Poales", "Caryophyllales", "Asterales", "Spermatophyta"))
newOrder <- master.ses.alpha[master.ses.alpha$clade == "Spermatophyta" & master.ses.alpha$metric == "mntd",] %>% arrange(ntaxa)
newOrder <- newOrder$summits
master.ses.alpha$summits <- factor(master.ses.alpha$summits, levels = newOrder)
master.ses.alpha$metric <- factor(master.ses.alpha$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

plot.pools <- ggplot(master.ses.alpha, aes(y=clade, x=summits, fill=obs.z))
plot.pools <- plot.pools + geom_tile(colour="white", alpha=.75)
plot.pools <- plot.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5, 5))
plot.pools <- plot.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.pools <- plot.pools + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.pools <- plot.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
#plot.pools <- plot.pools + facet_grid(pool ~ metric)#,space="free",scales="free", as.table = F)   
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
  axis.text.y = element_text(size = 10 * 0.8, colour = "black"),
  axis.text.x = element_text(size = 10 * 0.8, angle = -45, hjust = 0, colour = "black"),
  axis.text.x = element_text(size = 6 * 0.8, angle = -45, hjust = 0, colour = "black"),
  strip.text.x = element_text(size = 8, face="bold"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="figs/phylo_div_alpha_SES_tile.pdf")
plot.pools
dev.off()


############## Just species in talus (not meadow)
master.ses.alpha.talus <- read.csv(file="output/08_PhyloDiversity/MiSeq/alpha/static/alpine.phylogeny.pool.SES.talus.csv", row.names=1)
head(master.ses.alpha.talus)

master.ses.alpha.talus$clade <- factor(master.ses.alpha.talus$clade, levels = c("Saxifragales", "Poales", "Brassicales", "Caryophyllales", "Asterales", "Spermatophyta"))
newOrder <- master.ses.alpha.talus[master.ses.alpha.talus$clade == "Spermatophyta" & master.ses.alpha.talus$metric == "mntd",] %>% arrange(ntaxa)
newOrder <- newOrder$summits
master.ses.alpha.talus$summits <- factor(master.ses.alpha.talus$summits, levels = newOrder)
master.ses.alpha.talus$metric <- factor(master.ses.alpha.talus$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

plot.pools <- ggplot(master.ses.alpha.talus, aes(y=clade, x=summits, fill=obs.z))
plot.pools <- plot.pools + geom_tile(colour="white", alpha=.75)
plot.pools <- plot.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5, 5))
plot.pools <- plot.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.pools <- plot.pools + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.pools <- plot.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
#plot.pools <- plot.pools + facet_grid(pool ~ metric)#,space="free",scales="free", as.table = F)   
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
  axis.text.y = element_text(size = 10 * 0.8, colour = "black"),
  axis.text.x = element_text(size = 10 * 0.8, angle = -45, hjust = 0, colour = "black"),
  strip.text.x = element_text(size = 8, face="bold"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="figs/phylo_div_alpha_SES_tile.talus.pdf")
plot.pools
dev.off()


