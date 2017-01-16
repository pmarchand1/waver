
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
#' using the \code{\link[geosphere]{distGeo}} function from the geosphere R
#' package. All distance are expressed in meters.
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
#' @examples
#'  pt <- SpatialPoints(matrix(c(0, 0), ncol = 2),
#'                      proj4string = CRS("+proj=longlat"))
#'  # Shoreline is a rectangle from (-0.2, 0.25) to (0.3, 0.5)
#'  rect <- Polygon(cbind(c(rep(-0.2, 2), rep(0.3, 2), -0.2),
#'                        c(0.25, rep(0.3, 2), rep(0.25, 2))))
#'  land <- SpatialPolygons(list(Polygons(list(rect), ID = 1)),
#'                          proj4string = CRS("+proj=longlat"))
#'  fetch_len(pt, bearings = c(0, 45, 225, 315), land,
#'            dmax = 50000, spread = c(-10, 0, 10))
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
        in_water <- is.null(rgeos::gIntersects(p, shoreline, byid = TRUE,
                                               returnDense = FALSE)[[1]])
        if(!in_water) {
            warning("point on land, returning NA")
            return(setNames(rep(NA, length(bearings)), bearings))
        }
    }

    # Clip shoreline layer to a rectangle around point
    # to guarantee at least dmax on each side
    clip_rect <- get_clip_rect(p, dmax, projected)
    shore_clip <- tryCatch(
        rgeos::gIntersection(shoreline, clip_rect, byid = TRUE),
        # If it fails, try byid = FALSE
        error = function(e) tryCatch(
            rgeos::gIntersection(shoreline, clip_rect, byid = FALSE),
            error = function(e) {
                warning("Error clipping shoreline, returning NA")
                return(setNames(rep(NA, length(bearings)), bearings))
            }
        )
    )

    # Convert any polygons to lines to get line-line intersections later
    shore_clip <- convert_to_lines(shore_clip)
    # If no land within rectangle, return dmax for all bearings
    if (is.null(shore_clip)) {
        return(setNames(rep(dmax, length(bearings)), bearings))
    }

    # Calculate fetch
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
#' \code{fetch_len_multi} provides two methods to efficiently compute fetch length
#' for multiple points.
#'
#' With \code{method = "btree"}, the \code{\link[rgeos]{gBinarySTRtreeQuery}}
#' function from the rgeos package is called to determine which polygons in
#' \code{shoreline} could be within \code{dmax} of each point. This is a fast
#' calculation based on bounding box overlap.
#'
#' With \code{method = "clip"}, the \code{shoreline} layer is clipped to a polygon
#' formed by the union of rectangular buffers around each point.
#'
#' In both cases, \code{\link{fetch_len}} is then applied to each point,
#' using only the necessary portion of the shoreline.
#'
#' Generally, the "clip" method will produce the biggest time savings when
#' points are clustered within distances less than \code{dmax} (so their
#' clipping rectangles overlap), whereas the "btree" method will be more
#' efficient when the shoreline is composed of multiple polygons and points are
#' distant from each other.
#'
#' @param pts A SpatialPoints* object.
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline SpatialLines* or SpatialPolygons* object representing the
#'  shoreline.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing.
#' @param method Whether to use the "btree" (default) or "clip" method.
#'  See below for more details.
#' @param projected Should projected coordinates be used to calculate fetch?
#' @return A matrix of fetch lengths, with one row by point in \code{pts} and
#'  one column by bearing in \code{bearings}.
#' @seealso \code{\link{fetch_len}} for details on the fetch length computation.
#' @export
fetch_len_multi <- function(pts, bearings, shoreline, dmax,
                     spread = 0, method = c("btree", "clip"), projected = FALSE) {
    # Check inputs
    match.arg(method)
    if (!is(pts, "SpatialPoints")) stop("pts must be a SpatialPoints* object.")
    pts <- as(pts, "SpatialPoints")  # remove DataFrame part if there is one
    if (is(shoreline, "SpatialLines")) {
        shoreline <- as(shoreline, "SpatialLines")
    } else if (is(shoreline, "SpatialPolygons")) {
        shoreline <- as(shoreline, "SpatialPolygons")
    } else {
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

    # Create rectangular buffers around each point
    rect_list <- lapply(1:length(pts),
                        function(i) get_clip_rect(pts[i], dmax, projected))
    rect_buf <- do.call(rbind, c(rect_list, makeUniqueIDs = TRUE))

    if (method == "btree") {
        # Generate list of shoreline polygon IDs with bounding box overlap for each rectangle
        btree <- rgeos::gBinarySTRtreeQuery(shoreline, rect_buf)
        # Calculate fetch for point at index i using btree
        fetch_i <- function(i) {
            if (is.null(btree[[i]])) {
                setNames(rep(dmax, length(bearings)), bearings)
            } else {
                fetch_len(pts[i], bearings, shoreline[btree[[i]]], dmax,
                          spread, projected, check_inputs = FALSE)
            }
        }
        # Calculate fetch for all points and return a (points x bearings) matrix
        fetch_res <- t(vapply(1:length(pts), fetch_i, rep(0, length(bearings))))
    } else { # method == "clip"
        # Clip shoreline to a merged buffer around all points
        rect_buf <- rgeos::gUnaryUnion(rect_buf)
        sub_shore <- rgeos::gIntersection(shoreline, rect_buf, byid = TRUE)
        fetch_res <- t(
            vapply(1:length(pts),
                   function(i) fetch_len(pts[i], bearings, sub_shore, dmax,
                                         spread, projected, check_inputs = FALSE),
                   rep(0, length(bearings)))
        )
    }
    fetch_res
}


#### Helper functions below are not exported by the package ####

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
        min(geosphere::distGeo(p, inters))
    } else if (class(inters) == "SpatialLines") {
        min(geosphere::distGeo(p, lines_to_endpts(inters)))
    } else if (class(inters) == "SpatialCollections") {
        coord_mat <- rbind(coordinates(inters@pointobj),
                           coordinates(lines_to_endpts(inters@lineobj)))
        min(geosphere::distGeo(p, coord_mat))
    } else {
        warning(paste("Point at", c(p[1], p[2]),
                      "cannot calculate distance to shore, returning NA."))
        NA
    }
}

