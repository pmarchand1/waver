# File to test waver package

library(waver)

# Read shoreline layer (GSHSS_clip) and point of interest layer (Eco)
GSHHS <- rgdal::readOGR("/nfs/dgill-data/GIS/Contextual data/GSHHS_Shoreline/f", "GSHHS_f_L1")
GSHHS <- as(GSHHS, "SpatialPolygons")
GSHHS_clip <- rgdal::readOGR("/nfs/dgill-data/GIS/Contextual data/GSHHS_Shoreline", "Shoreline_50km_EcoClip_27May15")
GSHHS_clip <- as(GSHHS_clip, 'SpatialPolygons')
Eco <- rgdal::readOGR("/nfs/dgill-data/GIS/Ecological data", "Combined_Eco_data_28May15_landadj")
Eco <- coordinates(Eco)

# Bearings for which to calculate fetch
bearings <- seq(0, 337.5, 22.5)
spread <- seq(-10, 10, 2.5)  # relative sub-bearings to average over
#can also do every 22.5/16 degrees
# Maximum distance
dmax <- 50000

# Try calculating fetch for 1st point
system.time(fetch_res <- fetch_len(Eco[1, ], bearings, GSHHS_clip, dmax, spread))

# Multiple points (apply here vs. fetch_multi)
system.time(fetch10 <- t(apply(Eco[1:10, ], 1, fetch_len, bearings, GSHHS_clip, dmax, spread)))

system.time(fetchmult <- fetch_len_multi(Eco[1:10, ], bearings, GSHHS_clip, dmax, spread))


