#####################################################################################################################
############# Phylogenetic diversity within alpine summits (alpha) ################################################## 
############# Static Null Model ##################################################################################### 
############# Hannah E. Marx, 20 Jan 2019 ########################################################################### 
#####################################################################################################################

###### Statistical test of PD with elevation:

master.ses.alpha <- read.csv(file="output/10_PhyloDiversity/alpha/sawtooth.alpha.RD.SES.csv", row.names=1)
head(master.ses.alpha)
str(master.ses.alpha)
master.ses.alpha <- full_join(master.ses.alpha, cbind("summits" = rownames(sawMeta), sawMeta["Elevation"], sawMeta["Treeline..m."], sawMeta["Extent.elevation"], sawMeta["Range"]))
master.ses.alpha$sig
master.ses.alpha$metric <- factor(master.ses.alpha$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))
master.ses.alpha$data <- factor(master.ses.alpha$data, labels = c("total"="Total Data"))
#master.ses.alpha$clade <- factor(master.ses.alpha$clade, levels = c("Brassicales", "Ericales", "Caryophyllales", "Poales", "Asterales", "Tracheophyta"))
master.ses.alpha$clade <- factor(master.ses.alpha$clade, levels = c("Tracheophyta", "Asterales", "Poales", "Caryophyllales", "Lamiales", "Brassicales", "Ericales"))
master.ses.alpha$pool <- factor(master.ses.alpha$pool, levels = c("All Alpine", "Talus", "Meadow"))
master.ses.alpha$summits <- factor(master.ses.alpha$summits, levels = c("Horstmann Peak",
                                                                        "Braxon Peak",
                                                                        "Thompson Peak",
                                                                        "Snowyside Peak",
                                                                        "Mount Cramer",
                                                                        "D.O. Lee Peak",
                                                                        "Salzburger Spitzl", 
                                                                        "Castle Peak",
                                                                        "Hyndman Peak"))
master.ses.alpha$summits <- factor(master.ses.alpha$summits, labels = c("Horstmann Peak" = "Horstmann Peak (3155 m)",
                                                                        "Braxon Peak" = "Braxon Peak (3156 m)",
                                                                        "Thompson Peak" =  "Thompson Peak (3203 m)",
                                                                        "Snowyside Peak" = "Snowyside Peak (3246 m)",
                                                                        "Mount Cramer" ="Mount Cramer (3266 m)",
                                                                        "D.O. Lee Peak" = "D.O. Lee Peak (3457 m)",
                                                                        "Salzburger Spitzl" = "Salzburger Spitzl (3536 m)", 
                                                                        "Castle Peak" = "Castle Peak (3601 m)",
                                                                        "Hyndman Peak"  ="Hyndman Peak (3660 m)"))
master.ses.alpha$data <- factor(master.ses.alpha$data, levels = c("Total Data"))



#################### Maximum Elevation #################### 


#### MNTD TOTAL DATA: #YVAR ~ XVAR
lm.mntd.trach <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Tracheophyta") %>%
     select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mntd.trach) 
par(mfrow=c(2,2)) 
##### Test assumptions of parametric model:
plot(lm.mntd.trach, which=1:4)      

lm.mntd.ast <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Asterales") %>%
             select(summits, obs.z, Elevation), formula = obs.z ~ Elevation) 
summary(lm.mntd.ast) 
par(mfrow=c(2,2)) 
plot(lm.mntd.ast, which=1:4)      

lm.mntd.poa <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Poales") %>%
             select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mntd.poa) 
par(mfrow=c(2,2)) 
plot(lm.mntd.poa, which=1:4)  

lm.mntd.car <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Caryophyllales") %>%
             select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mntd.car) 
par(mfrow=c(2,2)) 
plot(lm.mntd.car, which=1:4)  

lm.mntd.lam <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Lamiales") %>%
             select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mntd.lam) 
par(mfrow=c(2,2)) 
plot(lm.mntd.lam, which=1:4)  

lm.mntd.bra <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Brassicales") %>%
                    select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mntd.bra) 
par(mfrow=c(2,2)) 
plot(lm.mntd.bra, which=1:4)  

lm.mntd.eri <- lm(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Ericales") %>%
             select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mntd.eri) 
par(mfrow=c(2,2)) 
plot(lm.mntd.eri, which=1:4)  



#### MPD TOTAL DATA
lm.mpd.trach <- lm(data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Tracheophyta") %>%
                      select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mpd.trach) 
par(mfrow=c(2,2)) 
plot(lm.mpd.trach, which=1:4)      

lm.mpd.ast <- lm(data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Asterales") %>%
                    select(summits, obs.z, Elevation), formula = obs.z ~ Elevation) 
summary(lm.mpd.ast) 
par(mfrow=c(2,2)) 
plot(lm.mpd.ast, which=1:4)      

lm.mpd.poa <- lm(data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Poales") %>%
                    select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mpd.poa) 
par(mfrow=c(2,2)) 
plot(lm.mpd.poa, which=1:4)  

lm.mpd.car <- lm(data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Caryophyllales") %>%
                    select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mpd.car) 
par(mfrow=c(2,2)) 
plot(lm.mpd.car, which=1:4)  

lm.mpd.lam <- lm(data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Lamiales") %>%
                    select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mpd.lam) 
par(mfrow=c(2,2)) 
plot(lm.mpd.lam, which=1:4)  


lm.mpd.bra <- lm(data = na.omit(filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Brassicales") %>%
                                  select(summits, obs.z, Elevation)), formula = obs.z ~ Elevation)
summary(lm.mpd.bra) 
par(mfrow=c(2,2)) 
plot(lm.mpd.bra, which=1:4)  


lm.mpd.eri <- lm(data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Ericales") %>%
                    select(summits, obs.z, Elevation), formula = obs.z ~ Elevation)
summary(lm.mpd.eri) 
par(mfrow=c(2,2)) 
plot(lm.mpd.eri, which=1:4)  



pdf(file = "figs/elevation/lm_elevationMax.pdf", width=4, height=10)
ggplot(filter(master.ses.alpha, pool == "All Alpine" & data == "Total Data"), aes(x=(Elevation), y=obs.z)) +
  geom_point() + 
  ylim(-4,4) +
  stat_smooth(mapping = aes(Elevation), method ="lm")  + 
  facet_grid(clade~metric) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position="none", 
        axis.text.x = element_text(angle = -45, hjust = -.5)) +
  labs( y = "SES", x = "Maximum elevation (m)")
dev.off()


######## Linear Mixed Model: fixed (max elevation) and all random effects (range -- has been sampled from an infinite populaition)
# https://dynamicecology.wordpress.com/2015/11/04/is-it-a-fixed-or-random-effect/
# https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/

# Fit a model including fixed (max elevation) and all random effects (range)
lmm.mntd.trac <- lmer(obs.z ~ Elevation + (1 | Range), data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Tracheophyta") %>%
             select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mntd.trac) 
#R2m        R2c
#[1,] 0.07837618 0.07837618
summary(lmm.mntd.trac)
plot(lmm.mntd.trac)

lmm.mntd.trac2 <- lme(data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Tracheophyta") %>%
      select(summits, obs.z, Elevation, Range), obs.z ~ Elevation, random = ~1|Range/Elevation)
coef(lmm.mntd.trac2) # quite low variation of coeff. with random effects
plot(random.effects(lmm.mntd.trac2)) # there should be symmetrical scatter around 0
plot(lmm.mntd.trac2) # much more heteroskedasticity 

# The marginal R squared values are those associated with your fixed effects, the conditional ones are those of your 
#fixed effects plus the random effects. Usually we will be interested in the marginal effects.

lmm.mntd.ast <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Asterales") %>%
                        select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mntd.ast) 
#0.01603439 0.02264945

lmm.mntd.car <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Caryophyllales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mntd.car) 
#0.129531 0.1338998

lmm.mntd.poa <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Poales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mntd.poa) 
# 0.001379589 0.006149563

lmm.mntd.lam <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Lamiales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mntd.lam) 
## NA

lmm.mntd.bra <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MNTD" & pool == "All Alpine" & data == "Total Data" & clade == "Brassicales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mntd.bra)
# 0.04582134 0.4426131


## MPD
lmm.mpd.trac <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Tracheophyta") %>%
                        select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mpd.trac) 
#R2m        R2c
# 0.4101094 0.4143077

lmm.mpd.ast <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Asterales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mpd.ast) 
# 0.001391309 0.6898869

lmm.mpd.car <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Caryophyllales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mpd.car) 
#0.218952 0.2222388

lmm.mpd.poa <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Poales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mpd.poa) 
# 0.2938108 0.2977528

lmm.mpd.lam <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Lamiales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mpd.lam) 
## NA

lmm.mpd.bra <- lmer(obs.z ~ Elevation + (Elevation | Range), data = filter(master.ses.alpha, metric == "MPD" & pool == "All Alpine" & data == "Total Data" & clade == "Brassicales") %>%
                       select(summits, obs.z, Elevation, Range))
r.squaredGLMM(lmm.mpd.bra)
# 5.86372e-05 0.003750683




########### Explainatory Distance Matrices : space and enviornmental (Elevation, Treeline..m.) distance 

spatialDist <- vegdist(sawMeta[c("WGS.N", "WGS.W")], method = "euclid")
enviroDistMaxElev <- vegdist(sawMeta["Elevation"], method = "euclid")
enviroDistRangeElev <- vegdist(sawMeta["Treeline..m."], method = "euclid")

assoSpace <- protest(X = decompo_beta$betadiv$SES_UniFrac_turn, Y = spatialDist, permutations = 10000)
summary(assoSpace)
plot(assoSpace)
assoSpace 

assoEnviroMaxElev <- protest(X = decompo_beta$betadiv$SES_UniFrac_turn, Y = enviroDistMaxElev, permutations = 10000)
plot(assoEnviroMaxElev)
assoEnviroMaxElev


## True turnover
MRM_spatialdist_turn <- MRM(decompo_beta$betadiv$SES_UniFrac_turn ~ spatialDist, nperm = 10000)
#$coef
#decompo_beta$betadiv$SES_UniFrac_turn   pval
#Int                                     1.3830650 0.2284
#spatialDist                            -0.5273684 0.4190

#$r.squared
#R2       pval 
#0.02275857 0.41900000 

#$F.test
#F    F.pval 
#0.7918118 0.4190000


MRM_envirodistMaxElev_turn <- MRM(decompo_beta$betadiv$SES_UniFrac_turn ~ enviroDistMaxElev, nperm = 10000)
#$coef
#decompo_beta$betadiv$SES_UniFrac_turn   pval
#Int                                        1.3484389729 0.2394
#enviroDistMaxElev                         -0.0008728403 0.4497

#$r.squared
#R2       pval 
#0.01470455 0.44970000 

#$F.test
#F    F.pval 
#0.5074159 0.4497000 

## PD component
MRM_spatialdist_pd <- MRM(decompo_beta$betadiv$SES_UniFrac_PD ~ spatialDist, nperm = 10000)
MRM_envirodistMaxElev_pd <- MRM(decompo_beta$betadiv$SES_UniFrac_PD ~ enviroDistMaxElev, nperm = 10000)

## Undecomposed Turnover
MRM_spatialdist <- MRM(decompo_beta$betadiv$SES_UniFrac ~ spatialDist, nperm = 10000)
MRM_envirodistMaxElev <- MRM(decompo_beta$betadiv$SES_UniFrac ~ enviroDistMaxElev, nperm = 10000)



