#####################################################################################################################
############# Statistical test of affect of alpine meadows on PD ################################################## 
############# Hannah E. Marx, 20 Jan 2019 ########################################################################### 
#####################################################################################################################

##### Filter data for summits with meadows
master.ses.alpha <- read.csv(file="output/10_PhyloDiversity/alpha/sawtooth.alpha.RD.SES.csv", row.names=1)

master.ses.alpha.meadows <- arrange(filter(master.ses.alpha, summits %in% c("Thompson Peak",
                                                                  "D.O. Lee Peak", 
                                                                  "Salzburger Spitzl",
                                                                  "Hyndman Peak")), clade)
master.ses.alpha.meadows$clade <- factor(master.ses.alpha.meadows$clade, levels = c( "Tracheophyta", "Asterales", "Poales","Caryophyllales",  "Lamiales", "Brassicales", "Ericales"))

##### Test assumptions of parametric models:
qqnorm(filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total" & pool == "Talus")$obs.z)
qqline(filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total" & pool == "Talus")$obs.z)

qqnorm(filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total" & pool == "All Alpine")$obs.z)
qqline(filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total" & pool == "All Alpine")$obs.z)

qqnorm(filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total" & pool == "Talus")$obs.z)
qqline(filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total" & pool == "Talus")$obs.z)

qqnorm(filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total" & pool == "All Alpine")$obs.z)
qqline(filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total" & pool == "All Alpine")$obs.z)

# Sample Variance: p-value > 0.05 = variences are homogeneous
var.test(filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total" & pool == "Talus")$obs.z, 
         filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total" & pool == "All Alpine")$obs.z)
#p-value = 0.7106

var.test(filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total" & pool == "Talus")$obs.z, 
         filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total" & pool == "All Alpine")$obs.z)
#p-value = 0.4213




pdf(file = "figs/habitat/box_habitat_MNTD.pdf", width=10)
ggplot(data = filter(master.ses.alpha.meadows, metric == "mntd" & pool != "Meadow" & data == "total"), 
       aes(x = pool, y = obs.z, fill = pool)) + geom_boxplot() + facet_grid(~clade) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs( y = "SES MNTD") +
  ylim(-2,2) 
dev.off()

pdf(file = "figs/habitat/box_habitat_MPD.pdf", width=10)
ggplot(data = filter(master.ses.alpha.meadows, metric == "mpd" & pool != "Meadow" & data == "total"), 
       aes(x = pool, y = obs.z, fill = pool)) + geom_boxplot() + facet_grid(~clade) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs( y = "SES MPD") +
  ylim(-2,2) 
dev.off()

#### MNTD
t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Tracheophyta") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Tracheophyta") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Asterales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Asterales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Poales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Poales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Caryophyllales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Caryophyllales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Lamiales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Lamiales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Brassicales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Brassicales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mntd" & pool == "All Alpine" & data == "total" & clade == "Ericales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mntd" & pool == "Talus" & data == "total" & clade == "Ericales") %>% select(obs.z))


#### mpd
t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Tracheophyta") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Tracheophyta") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Asterales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Asterales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Poales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Poales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Caryophyllales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Caryophyllales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Lamiales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Lamiales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Brassicales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Brassicales") %>% select(obs.z))

t.test(x = filter(master.ses.alpha.meadows, metric == "mpd" & pool == "All Alpine" & data == "total" & clade == "Ericales") %>% select(obs.z), 
       filter(master.ses.alpha.meadows, metric == "mpd" & pool == "Talus" & data == "total" & clade == "Ericales") %>% select(obs.z))


