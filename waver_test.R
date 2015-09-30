# File to test waver package

library(fetch)

# Read shoreline layer (GSHSS_clip) and point of interest layer (Eco)
# GSHHS <- readOGR("/nfs/dgill-data/GIS/Contextual data/GSHHS_Shoreline/f", "GSHHS_f_L1")
GSHHS_clip <- rgdal::readOGR("/nfs/dgill-data/GIS/Contextual data/GSHHS_Shoreline", "Shoreline_50km_EcoClip_27May15")
Eco <- readOGR("/nfs/dgill-data/GIS/Ecological data", "Combined_Eco_data_28May15_landadj")


# Bearings for which to calculate fetch
bearings <- seq(0, 337.5, 22.5)
spread <- seq(-10, 10, 2.5)  # relative sub-bearings to average over
#can also do every 22.5/16 degrees
# Maximum distance
dmax <- 50000

coords <- data.frame(coordinates(Eco))
colnames(coords) <- c('long', 'lat')

GSHHS_clip <- as(GSHHS_clip, 'SpatialPolygons')
