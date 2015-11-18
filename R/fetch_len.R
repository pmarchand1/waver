
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
#' If \code{projected} is FALSE (the default), the input data must be in WGS84
#' geographic (longitude, latitude) coordinates. Geodesic distances are calculated
#' using the Vincenty ellipsoid formula from the geosphere R package. All
#' distance are expressed in meters.
#'
#' If \code{projected} is TRUE, the input data (\code{p} and \code{shoreline})
#' must share the same projection. Projected distances are calculated with the
#' rgeos R package. All distances are expressed in the projection's coordinates.
#'
#' If the shoreline layer is given as SpatialPolygons*, the function verifies
#' that the input point is outside all polygons (i.e. in water). If this is
#' not the case, it issues a warning and returns a vector of \code{NA}.
#'
#' @param p SpatialPoints* object of length 1 (single point).
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline SpatialLines* or SpatialPolygons* object representing the
#'  shoreline.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing (see details).
#' @param projected Should projected coordinates be used to calculate fetch?
#' @param check_inputs Should the validity of inputs be checked? It is
#' recommended to keep this TRUE, unless this function is called repeatedly from
#' another function that already checks inputs.
#' @return A named vector representing the fetch length for each direction
#'  given in \code{bearings}.
#' @seealso \code{\link{fetch_len_multi}} for an efficient alternative when
#'  computing fetch length for multiple points.
#' @export
fetch_len <- function(p, bearings, shoreline, dmax,
                      spread = 0, projected = FALSE, check_inputs = TRUE) {
    if (check_inputs) {
        if (!is(p, "SpatialPoints")) stop("p must be a SpatialPoints* object.")
        p <- as(p, "SpatialPoints")  # remove DataFrame part if there is one
        if(length(p) != 1) stop("p must be a single point.")
        if (!(is(shoreline, "SpatialLines") || is(shoreline, "SpatialPolygons"))) {
            stop("shoreline must be a SpatialLines* or SpatialPolygons* object.")
        }
        if (projected) {
            if (!is.projected(p) || !is.projected(shoreline)) {
                stop("cannot use long/lat coordinates if projected = TRUE.")
            }
            if (proj4string(p) != proj4string(shoreline)) {
                stop("projections of p and shoreline do not match.")
            }
        } else if (is.projected(p) || is.projected(shoreline)) {
                stop(paste("p and shoreline must have unprojected (long/lat)",
                           "coordinates if projected = FALSE."))
        }
        if (!is.vector(bearings, "numeric")) stop("bearings must be a numeric vector.")
        if (!is.vector(spread, "numeric")) stop("spread must be a numeric vector.")
        if (!is.vector(dmax, "numeric") || length(dmax) != 1 || dmax <= 0) {
            stop("dmax must be a single number greater than 0.")
        }
    }

    # If shoreline is a polygons (land) layer, check that point is not on land
    if (is(shoreline, "SpatialPolygons")) {
        # Note: 'as' only to remove DataFrame part of Spatial objects
        in_water <- is.na(over(as(p, "SpatialPoints"),
                               as(shoreline, "SpatialPolygons")))
        if(!in_water) {
            warning("point on land, returning NA")
            return(setNames(rep(NA, length(bearings)), bearings))
        }
    }

    # Clip shoreline layer to a rectangle around point
    # to guarantee at least dmax on each side
    clip_rect <- get_clip_rect(p, dmax, projected)
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
        fetch_res <- vapply(bearings,
                   function(b) dist_shore(p, shore_clip, b, dmax, projected), 0)
    } else {
        # calculate the distance to shore for each sub-bearing
        bear_mat <- outer(bearings, spread, "+")
        dists <- vapply(bear_mat,
                   function(b) dist_shore(p, shore_clip, b, dmax, projected), 0)
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


#' Calculate the fetch length for multiple points
#'
#' \code{fetch_len_multi} provides an efficient method to compute fetch length
#' for multiple points.
#'
#' This function clips the \code{shoreline} layer to a rectangle around each
#' point in \code{pts} before applying \code{\link{fetch_len}} to all points.
#' This saves computation time if points of interest are spatially clustered
#' and their rectangular buffers overlap.
#'
#' @param pts A SpatialPoints* object.
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline SpatialLines* or SpatialPolygons* object representing the
#'  shoreline.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing.
#' @param projected Should projected coordinates be used to calculate fetch?
#' @return A matrix of fetch lengths, with one row by point in \code{pts} and
#'  one column by bearing in \code{bearings}.
#' @seealso \code{\link{fetch_len}} for details on the fetch length computation.
#' @export
fetch_len_multi <- function(pts, bearings, shoreline, dmax,
                            spread = 0, projected = FALSE) {
    # Check inputs
    if (!is(pts, "SpatialPoints")) stop("pts must be a SpatialPoints* object.")
    pts <- as(pts, "SpatialPoints")  # remove DataFrame part if there is one
    if (!(is(shoreline, "SpatialLines") || is(shoreline, "SpatialPolygons"))) {
        stop("shoreline must be a SpatialLines* or SpatialPolygons* object.")
    }
    if (projected) {
        if (!is.projected(pts) || !is.projected(shoreline)) {
            stop("cannot use long/lat coordinates if projected = TRUE.")
        }
        if (proj4string(pts) != proj4string(shoreline)) {
            stop("projections of pts and shoreline do not match.")
        }
    } else if (is.projected(pts) || is.projected(shoreline)) {
            stop(paste("pts and shoreline must have unprojected (long/lat)",
                       "coordinates if projected = FALSE."))
    }
    if (!is.vector(bearings, "numeric")) stop("bearings must be a numeric vector.")
    if (!is.vector(spread, "numeric")) stop("spread must be a numeric vector.")
    if (!is.vector(dmax, "numeric") || length(dmax) != 1 || dmax <= 0) {
        stop("dmax must be a single number greater than 0.")
    }

    # Clip shoreline to a merged buffer around all points
    rect_list <- lapply(1:length(pts),
                        function(i) get_clip_rect(pts[i], dmax, projected))
    rect_buf <- do.call(rbind, c(rect_list, makeUniqueIDs = TRUE))
    rect_buf <- rgeos::gUnaryUnion(rect_buf)
    sub_shore <- rgeos::gIntersection(shoreline, rect_buf, byid = TRUE)

    # Calculate fetch for all points and return a (points x bearings) matrix
    fetch_res <- t(
        vapply(1:length(pts),
               function(i) fetch_len(pts[i], bearings, shoreline, dmax,
                                     spread, projected, check_inputs = FALSE),
               rep(0, length(bearings)))
        )
    fetch_res
}


#### Helper functions below are not exported by the package ####


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
            clip_rect <- SpatialPolygons(list(Polygons(list(westr), ID = 1),
                                              Polygons(list(eastr), ID = 2)),
                                         proj4string = CRS(proj4string(p)))
        } else if(long + xbuf > 180) {
            westr <- poly_rect(-180, lat - ybuf, long + xbuf - 360, lat + ybuf)
            eastr <- poly_rect(long - xbuf, lat - ybuf, 180, lat + ybuf)
            clip_rect <- SpatialPolygons(list(Polygons(list(westr), ID = 1),
                                              Polygons(list(eastr), ID = 2)),
                                         proj4string = CRS(proj4string(p)))
        } else {
            r1 <- poly_rect(long - xbuf, lat - ybuf, long + xbuf, lat + ybuf)
            clip_rect <- SpatialPolygons(list(Polygons(list(r1), ID = 1)),
                                         proj4string = CRS(proj4string(p)))
        }
    }
    clip_rect
}

# Returns the distance from point p to shoreline
# following bearing bear and up to distance dmax
# all distances in projection units, or meters if projected = FALSE
dist_shore <- function(p, shoreline, bear, dmax, projected) {
    if (projected) {
        # Draw line of length dmax with given start point and bearing
        bline <- bearing_line(p, bear, dmax)
        # Return (minimum) distance from p1 to intersection of geo_line and shoreline
        # If no intersection, fetch is dmax
        land_int <- rgeos::gIntersection(bline, shoreline)
        if (is.null(land_int)) {
            dmax
        } else {
            rgeos::gDistance(p, land_int)
        }
    } else {
        # Draw geodesic line of length dmax with given start point and bearing
        # Line drawn with a point every dmax/500
        geo_line <- geosphere::gcIntermediate(p, geosphere::destPoint(p, bear, dmax),
                     n = 500, sp = TRUE, breakAtDateLine = TRUE, addStartEnd = TRUE)
        geo_line <- spTransform(geo_line, CRS(proj4string(shoreline)))
        # Return (minimum) distance from p to intersection of geo_line and shoreline
        # If no intersection, fetch is dmax
        land_int <- rgeos::gIntersection(geo_line, shoreline)
        if (is.null(land_int)) {
            dmax
        } else {
            dist_min(p, land_int)
        }
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

