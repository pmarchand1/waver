#### Some utility functions to create / manipulate spatial objects (not exported)

# Create a Polygon object corresponding to rectangle with given coords
poly_rect <- function(xmin, ymin, xmax, ymax) {
    Polygon(cbind(c(rep(xmin, 2), rep(xmax, 2), xmin),
                  c(ymin, rep(ymax, 2), rep(ymin, 2))))
}

# SpatialLines to endpoint coordinates
lines_to_endpts <- function(slines) {
    endpts_list <- lapply(coordinates(slines),
                          function(x) matrix(unlist(x), 2, 2))
    do.call(rbind, endpts_list)
}

# Create a SpatialLines starting at p, of length len,
#  along given bearing (in degrees)
# Note: This function is meant to be used with projected coordinates
bearing_line <- function(p, bearing, len) {
    l1 <- Line(rbind(
        coordinates(p),
        coordinates(p) + len * c(sinpi(bearing / 180), cospi(bearing / 180))
    ))
    SpatialLines(list(Lines(list(l1), ID = 1)), CRS(proj4string(p)))
}

# If sp_obj is a SpatialPolygons or SpatialCollections, convert to SpatialLines
# If it has no polygons or lines (e.g. SpatialPoints), return NULL
convert_to_lines <- function(sp_obj) {
    res <- NULL
    if (is(sp_obj, "SpatialLines")) {
        res <- sp_obj
    } else if (is(sp_obj, "SpatialPolygons")) {
        res <- as(sp_obj, "SpatialLines")
    } else if (is(sp_obj, "SpatialCollections")) {
        if (!is.null(sp_obj@polyobj) && !is.null(sp_obj@lineobj)) {
            res <- rbind(as(sp_obj@polyobj, "SpatialLines"), sp_obj@lineobj,
                         makeUniqueIDs = TRUE)
        } else if (!is.null(sp_obj@polyobj)) {
            res <- as(sp_obj@polyobj, "SpatialLines")
        } else if (!is.null(sp_obj@lineobj)) {
            res <- sp_obj@lineobj
        }
    }
    res
}

# Create clipping rectangle around point p (SpatialPoints of length 1)
#  to guarantee at least dmax on each side
# dmax either in meters (if projected = FALSE) or in the projection's coordinates
get_clip_rect <- function(p, dmax, projected) {
    if (projected) {
        x <- coordinates(p)[1]
        y <- coordinates(p)[2]
        r1 <- poly_rect(x - dmax, y - dmax, x + dmax, y + dmax)
        clip_rect <- SpatialPolygons(list(Polygons(list(r1), ID = 1)),
                                     proj4string = CRS(proj4string(p)))
    } else {
        lat_dist <- 111600 # approx. distance (in m) between degrees of latitude
        long <- coordinates(p)[1]
        lat <- coordinates(p)[2]
        ybuf = dmax / lat_dist
        xbuf = ybuf / cospi(abs(lat) / 180)
        # Split clip_rect in two if it would overlap international date line
        if (long - xbuf < -180) {
            westr <- poly_rect(-180, lat - ybuf, long + xbuf, lat + ybuf)
            eastr <- poly_rect(long - xbuf + 360, lat - ybuf, 180, lat + ybuf)
            clip_rect <- SpatialPolygons(list(Polygons(list(westr, eastr), ID = 1)),
                                         proj4string = CRS(proj4string(p)))
        } else if(long + xbuf > 180) {
            westr <- poly_rect(-180, lat - ybuf, long + xbuf - 360, lat + ybuf)
            eastr <- poly_rect(long - xbuf, lat - ybuf, 180, lat + ybuf)
            clip_rect <- SpatialPolygons(list(Polygons(list(westr, eastr), ID = 1)),
                                         proj4string = CRS(proj4string(p)))
        } else {
            r1 <- poly_rect(long - xbuf, lat - ybuf, long + xbuf, lat + ybuf)
            clip_rect <- SpatialPolygons(list(Polygons(list(r1), ID = 1)),
                                         proj4string = CRS(proj4string(p)))
        }
    }
    clip_rect
}

