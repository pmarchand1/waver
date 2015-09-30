#library(rgdal) # to load file
#library(geosphere) # for great circle distance calculations
#library(rgeos) # for spatial intersection
library(sp)

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

plonglat <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#########################
# All the functions below
#########################

# Create clipping rectangle around point
# to guarantee at least dmax on each side
get_clip_rect <- function(long, lat, dmax) {
    lat_dist <- 111600 # approx. distance (in m) between degrees of latitude
    ybuf = dmax / lat_dist
    xbuf = ybuf / cospi(abs(lat) / 180)
    if (long - xbuf < -180) {
        westr <- poly_rect(-180, lat - ybuf, long + xbuf, lat + ybuf)
        eastr <- poly_rect(long - xbuf + 360, lat - ybuf, 180, lat + ybuf)
        clip_rect <- SpatialPolygons(list(Polygons(list(westr), ID = 1),
                                          Polygons(list(eastr), ID = 2)),
                                     proj4string = CRS(plonglat))
    } else if(long + xbuf > 180) {
        westr <- poly_rect(-180, lat - ybuf, long + xbuf - 360, lat + ybuf)
        eastr <- poly_rect(long - xbuf, lat - ybuf, 180, lat + ybuf)
        clip_rect <- SpatialPolygons(list(Polygons(list(westr), ID = 1),
                                          Polygons(list(eastr), ID = 2)),
                                     proj4string = CRS(plonglat))
    } else {
        r1 <- poly_rect(long - xbuf, lat - ybuf, long + xbuf, lat + ybuf)
        clip_rect <- SpatialPolygons(list(Polygons(list(r1), ID = 1)),
                                     proj4string = CRS(plonglat))
    }
    clip_rect
}


# This function calculates the fetch lengths (in meters)
# around point (long, lat) for a set of bearings, with maximum dmax
# Optionally, may include a vector spread representing sub-bearings to average
#  fetch over (sub-bearing directions are relative to main bearing)
fetch_len <- function(long, lat, bearings, shoreline, dmax, spread = 0) {
    # Clip shoreline layer to a rectangle around point
    # to guarantee at least dmax on each side
    clip_rect <- get_clip_rect(long, lat, dmax)
    land_clip <- rgeos::gIntersection(GSHHS_clip, clip_rect, byid = TRUE)
    # If no land within rectangle, return dmax for all bearings
    if (is.null(land_clip)) return(rep(dmax, length(bearings)))
    # Convert shoreline from polygon to line (to get line-line intersections below)
    land_clip <- as(land_clip, 'SpatialLines')

    if (all(spread == 0)) {
        # if no sub-bearings, just return distance to shore for each bearing
        sapply(bearings, function(b) dist_shore(c(long, lat), land_clip, b, dmax))
    } else {
        # calculate the distance to shore for each sub-bearing
        bear_mat <- outer(bearings, spread, "+")
        dists <- sapply(bear_mat,
                        function(b) dist_shore(c(long, lat), land_clip, b, dmax))
        dim(dists) <- dim(bear_mat)
        # return weighted means of the sub-bearing fetch values
        #  with weights proportional to the cosine (relative to their main bearing)
        weights <- cospi(spread / 180)
        weights <- weights / sum(weights)
        dists %*% weights
    }
}


# Returns the distance from point p to shoreline
# following bearing bear and up to distance dmax
dist_shore <- function(p, shoreline, bear, dmax) {
    # Draw geodesic line of length dmax with given start point and bearing
    # Line drawn with a point every dmax/500
    geo_line <- gcIntermediate(p, destPoint(p, bear, dmax), n = 500, sp = TRUE,
                               breakAtDateLine = TRUE, addStartEnd = TRUE)
    geo_line <- spTransform(geo_line, CRS(proj4string(GSHHS_clip)))
    # Return (minimum) distance from p1 to intersection of geo_line and shoreline
    # If no intersection, fetch is dmax
    land_int <- gIntersection(geo_line, shoreline)
    if (is.null(land_int)) {
        dmax
    } else {
        dist_min(p, land_int)
    }
}


# Returns the minimum distance between the focal point p and inters, the result
#  of an intersection between SpatialLines which could include points and lines
dist_min <- function(p, inters) {
    if (class(inters) == 'SpatialPoints') {
        min(distVincentyEllipsoid(p, inters))
    } else if (class(inters) == 'SpatialLines') {
        min(distVincentyEllipsoid(p, lines_to_endpts(inters)))
    } else if (class(inters) == 'SpatialCollections') {
        coord_mat <- rbind(coordinates(inters@pointobj),
                           coordinates(lines_to_endpts(inters@lineobj)))
        min(distVincentyEllipsoid(p, coord_mat))
    } else {
        warning(paste('Point at', long,',', lat,
                      'cannot calculate distance to shore, returning NA.'))
        NA
    }
}

