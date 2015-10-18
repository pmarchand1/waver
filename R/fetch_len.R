### TODO: check inputs to fetch_len, check if any point is in water


#' Calculate the fetch length around a point
#'
#' Given a point, a shoreline layer and a vector of wind directions (bearings),
#' \code{fetch_len} calculates the distance from point to shore for each bearing.
#'
#' The fetch length (or fetch) is the distance of open water over which the wind
#' can blow in a specific direction. Note that bearings represent the direction
#' from where the wind originates.
#'
#' The optional \code{spread} argument defines relative directions that are
#' added to each main bearing to produce a set of sub-bearings. The fetch lengths
#' calculated for each sub-bearing are averaged with weights proportional to
#' \code{cos(spread)}. By default, \code{spread = 0} and fetch length is
#' calculated for the main bearings only.
#'
#' This function converts all Spatial* objects to WGS84 geographic
#' (longitude, latitude) coordinates and calculates geodesic distances using
#' the Vincenty ellipsoid formula (from the geosphere R package).
#'
#' If the shoreline layer is given as SpatialPolygons*, the function verifies
#' that the input point is outside all polygons (i.e. in water). If this is
#' not the case, it issues a warning and returns a vector of \code{NA}.
#'
#' @param p Coordinates of the point from which to calculate fetch, in degrees.
#'  Either a vector of two numbers, a matrix of two columns (longitude,
#'  latitude), or a SpatialPoints* object.
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline SpatialLines* or SpatialPolygons* object representing the
#'  shoreline.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing (see details).
#' @return A named vector representing the fetch length for each direction
#'  given in \code{bearings}.
#'
#' @export
fetch_len <- function(p, bearings, shoreline, dmax, spread = 0) {
    # Convert p to vector if necessary
    if (is.matrix(p)) p <- as.vector(p)
    if (is(p, "SpatialPoints")) p <- as.vector(coordinates(p))
    # Clip shoreline layer to a rectangle around point
    # to guarantee at least dmax on each side
    clip_rect <- get_clip_rect(p, dmax)
    shore_clip <- rgeos::gIntersection(shoreline, clip_rect, byid = TRUE)
    # If no land within rectangle, return dmax for all bearings
    if (is.null(shore_clip) || is(shore_clip, "SpatialPoints")) {
        return(rep(dmax, length(bearings)))
    }

    # Convert any polygons to lines to get line-line intersections later
    if (is(shore_clip, "SpatialPolygons")) {
        shore_clip <- as(shore_clip, "SpatialLines")
    } else if (is(shore_clip, "SpatialCollections")) {
        if (!is.null(shore_clip@polyobj)) {
            shore_clip <- rbind(as(shore_clip@polyobj, "SpatialLines"),
                                shore_clip@lineobj)
        } else if (!is.null(shore_clip@lineobj)) {
            shore_clip <- shore_clip@lineobj
        } else {  # no land in buffer
            return(rep(dmax, length(bearings)))
        }
    }

    if (all(spread == 0)) {
        # if no sub-bearings, just return distance to shore for each bearing
        fetch_res <- sapply(bearings,
                            function(b) dist_shore(p, shore_clip, b, dmax))
    } else {
        # calculate the distance to shore for each sub-bearing
        bear_mat <- outer(bearings, spread, "+")
        dists <- sapply(bear_mat,
                        function(b) dist_shore(p, shore_clip, b, dmax))
        dim(dists) <- dim(bear_mat)
        # return weighted means of the sub-bearing fetch values
        #  with weights proportional to the cosine (relative to their main bearing)
        weights <- cospi(spread / 180)
        weights <- weights / sum(weights)
        fetch_res <- as.vector(dists %*% weights)
    }

    names(fetch_res) <- as.character(bearings)
    fetch_res
}


# Calculate the fetch length for multiple points
# Right now pts can only be a matrix
#' @export
fetch_len_multi <- function(pts, bearings, shoreline, dmax, spread = 0) {
    # Clip shoreline to a merged buffer around all points
    rect_list <- apply(pts, 1, get_clip_rect, dmax)
    rect_buf <- do.call(rbind, c(rect_list, makeUniqueIDs = TRUE))
    rect_buf <- rgeos::gUnaryUnion(rect_buf)
    sub_shore <- rgeos::gIntersection(shoreline, rect_buf, byid = TRUE)

    # Calculate fetch for all points and return a (points x bearings) matrix
    fetch_res <- t(apply(pts, 1, fetch_len, bearings, sub_shore, dmax, spread))
    colnames(fetch_res) <- as.character(bearings)
    fetch_res
}


# Helper functions below are not exported by the package

# Constant for longlat WGS projection
plonglat <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Create clipping rectangle around point p (long, lat)
#   to guarantee at least dmax (in m) on each side
get_clip_rect <- function(p, dmax) {
    lat_dist <- 111600 # approx. distance (in m) between degrees of latitude
    long <- p[1]
    lat <- p[2]
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


# Returns the distance from point p to shoreline
# following bearing bear and up to distance dmax
dist_shore <- function(p, shoreline, bear, dmax) {
    # Draw geodesic line of length dmax with given start point and bearing
    # Line drawn with a point every dmax/500
    geo_line <- geosphere::gcIntermediate(p, geosphere::destPoint(p, bear, dmax),
                 n = 500, sp = TRUE, breakAtDateLine = TRUE, addStartEnd = TRUE)
    geo_line <- spTransform(geo_line, CRS(proj4string(shoreline)))
    # Return (minimum) distance from p1 to intersection of geo_line and shoreline
    # If no intersection, fetch is dmax
    land_int <- rgeos::gIntersection(geo_line, shoreline)
    if (is.null(land_int)) {
        dmax
    } else {
        dist_min(p, land_int)
    }
}


# Returns the minimum distance between the focal point p and inters, the result
#  of an intersection between SpatialLines which could include points and lines
dist_min <- function(p, inters) {
    if (class(inters) == "SpatialPoints") {
        min(geosphere::distVincentyEllipsoid(p, inters))
    } else if (class(inters) == "SpatialLines") {
        min(geosphere::distVincentyEllipsoid(p, lines_to_endpts(inters)))
    } else if (class(inters) == "SpatialCollections") {
        coord_mat <- rbind(coordinates(inters@pointobj),
                           coordinates(lines_to_endpts(inters@lineobj)))
        min(geosphere::distVincentyEllipsoid(p, coord_mat))
    } else {
        warning(paste("Point at", c(p[1], p[2]),
                      "cannot calculate distance to shore, returning NA."))
        NA
    }
}

