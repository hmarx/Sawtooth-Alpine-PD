
source("R/addScaleBarMap.R")

## Map resources:
#https://lib.stanford.edu/gis-branner-library/gis-data-north-america
#http://pubs.usgs.gov/of/2005/1235/
#https://insideidaho.org/webapps/search/iso_browse.aspx

# Geology Layer
#http://pubs.usgs.gov/of/2005/1305/

###### Forest Boundary
### Found National Forest layer on: http://data.fs.usda.gov/geodata/edw/datasets.php?dsetCategory=boundaries
## edited in qgis to include only Sawtooth National Forest -> saved layer 
layerName <- "SawtoothNF"  
data_projected <- readOGR(dsn="output/09_Maps/", layer=layerName) 
plot(data_projected)# +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs
projection(data_projected) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"

###### Transform shapefile to dataframe ###### 
# add to data a new column termed "id" composed of the rownames of data
data_projected@data$id <- rownames(data_projected@data)

# reproject the data onto a "longlat" projection
subsetTransform <- spTransform(data_projected, CRS("+proj=longlat"))

# create a data.frame from our spatial object
sawtoothPoints <- fortify(subsetTransform, region = "id")
head(sawtoothPoints)

# merge the "fortified" data with the data from our spatial object
sawtoothDF <- merge(sawtoothPoints, subsetTransform@data, by = "id")
head(sawtoothDF)
unique(sawtoothDF$FORESTNAME)

################################## Map of US with Idaho highlighted 
##################### PLOT ##################### 
usa <- map_data("state")
head(usa)

ID <- subset(usa, region %in% "idaho")
id.map <- ggplot() + geom_polygon(data = ID, aes(x = long, y = lat, group = group)) + coord_map()
#P + scaleBar(lon = -124, lat = 45, distanceLon = 50, distanceLat = 10, distanceLegend = 20, dist.unit = "km", orientation = FALSE)
id.map <- id.map + scaleBar(lon = -113, lat = 47.5, distanceLon = 50, distanceLat = 10, distanceLegend = 20, dist.unit = "km", 
                            arrow.length = 50, arrow.distance = 30, arrow.North.size = 4)
id.map <- id.map + theme(panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
                         panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_blank(),
                         axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank(), axis.title = element_blank(),
                         rect = element_blank(),
                         plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines"))
pdf(file="figs/maps/id.general.pdf")
id.map
dev.off()

usa.map <- map_data("state")
usa.map <- ggplot() + geom_polygon(data = usa.map, aes(x = long, y = lat, group = group), fill="NA", color="grey") + coord_map()
usa.map <- usa.map + geom_polygon(data = ID, aes(x = long, y = lat, group = group), fill="grey") 
usa.map <- usa.map + scaleBar(lon = -130, lat = 26, distanceLon = 500, distanceLat = 100, distanceLegend = 200, dist.unit = "km")
usa.map <- usa.map + theme(panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
                           panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_blank(),
                           axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(), axis.title = element_blank(),
                           rect = element_blank(),
                           plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines"))
pdf(file="figs/maps//usa.general.pdf")
usa.map
dev.off()



##################   Map of Ecrins NP zoomed in
########### Import custom raster basemap (made in qgis)

########### Import custom raster basemap (made in qgis)
### Baselayer from:
#http://felix.rohrba.ch/en/2016/awesome-basemap-layer-for-your-qgis-project/
#http://www.naturalearthdata.com/features/

layer <- stack("output/09_Maps/SawtoothBase.tiff") #EPSG: 3857
projection(layer) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"

# Reproject layer to lat long
sr <-"+proj=longlat"
projected_raster <- projectRaster(layer, crs = sr)
plot(projected_raster$SawtoothBase.1, col= grey(1:99/100))

## bounding box
domain <- c(-115.5, -113.5, 43.2, 44.45)

## crop layer to bounding box
sawtooth.base.crop <- crop(projected_raster, y=extent(domain))
plot(sawtooth.base.crop$SawtoothBase.1, col= grey(1:99/100))

## create raster dataframe
rast.table <- data.frame(xyFromCell(sawtooth.base.crop, 1:ncell(sawtooth.base.crop)), getValues(sawtooth.base.crop/255))
head(rast.table)
  
ses.all <- read.csv("output/08_PhyloDiversity/MiSeq//alpha//static//alpine.phylogeny.pool.SES.csv")
head(ses.all)
ses.all.meta <- merge(ses.all, sawMeta, by.x=10, by.y=0)
head(ses.all.meta)
ses.all.meta

cols = c("Hyndman Peak" = "#F8766D",
            "D.O. Lee Peak" = "#E58700",
            "Salzburger Spitzl" = "#A3A500",
            "Castle Peak" = "#00BA38",
            "Thompson Peak"= "#00C0AF",
            "Mount Cramer" = "#00B0F6",
            "Snowyside Peak"= "#B983FF",
            "Horstmann Peak"= "#E76BF3",
            "Braxon Peak"= "#FF67A4")

## plot summits by elevation
ggplot(data = rast.table, aes(x = x, y = y)) +
  geom_raster(aes(fill=SawtoothBase.1)) +
  scale_fill_gradientn(colours=c("grey61","grey100"), guide = "none") +
  geom_polygon(data=sawtoothDF, aes(x=long, y=lat, group=group), colour="white", fill="grey10", alpha=0.4) + 
  geom_point(data =filter(ses.all.meta, clade=="Spermatophyta" & metric=="mntd"), aes(x=WGS.W, y=WGS.N, 
                        colour = summits, size = Elevation), pch=17) + #size=6,  
  scale_colour_manual(values = cols) +
  scale_alpha_discrete(range=c(1,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('long') + ylab('lat') + 
  #theme(legend.position = "none") +
  scaleBar(lon = -114, lat = 44.2, distanceLon = 10, distanceLat = 3, distanceLegend =6, dist.unit = "km", 
           arrow.length = 12, arrow.distance = 10, arrow.North.size = 6) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  coord_equal() 
ggsave("figs/maps/summits.elev.pdf", width=8, height = 8)



ses.all.meta$metric <- factor(ses.all.meta$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

## plot summits by SES all alpine
ggplot(data = rast.table, aes(x = x, y = y)) +
  geom_raster(aes(fill=SawtoothBase.1)) +
  scale_fill_gradientn(colours=c("grey61","grey100"), guide = "none") +
  geom_polygon(data=sawtoothDF, aes(x=long, y=lat, group=group), colour="white", fill="grey10", alpha=0.4) + 
  geom_point(data =filter(ses.all.meta, clade=="Spermatophyta"), aes(x=WGS.W, y=WGS.N, 
                colour = obs.z), pch=17, size=6) + # & metric=="mntd"
  scale_color_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5, 5), guides(title = "SES")) +
  geom_point(data =filter(ses.all.meta, clade=="Spermatophyta"), aes(x=WGS.W, y=WGS.N, 
                       size=as.factor(sig))) +  # & metric=="mntd"
  scale_size_manual(values=c("1"=1, "0"=NA), guide="none") +
  facet_grid(metric ~ .) +
  scale_alpha_discrete(range=c(1,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('long') + ylab('lat') + 
  theme(strip.text.x = element_text(size = 10, face="bold"),
          strip.text.y = element_text(size = 10, face="bold")) +
  #theme(legend.position = "none") +
  #scaleBar(lon = -114, lat = 44.2, distanceLon = 10, distanceLat = 3, distanceLegend =6, dist.unit = "km", 
  #         arrow.length = 12, arrow.distance = 10, arrow.North.size = 6) +
  #scale_size_continuous(range = c(3, 8), guide = "none") +
  coord_equal() 
ggsave("figs/maps/summits.SESalpine.pdf", width=8, height = 10)
 


