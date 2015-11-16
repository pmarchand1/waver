
#' Calculate the fetch length around a point (projected coordinates)
#'
#' Given a point, a shoreline layer and a vector of wind directions (bearings),
#' \code{fetch_len_proj} calculates the distance from point to shore for each bearing.
#' This version of the function assumes the point and shoreline share the same projection.
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
#' If the shoreline layer is given as SpatialPolygons*, the function verifies
#' that the input point is outside all polygons (i.e. in water). If this is
#' not the case, it issues a warning and returns a vector of \code{NA}.
#'
#' @param p SpatialPoints* object of length 1 (single point).
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline SpatialLines* or SpatialPolygons* object representing the
#'  shoreline, in the same projection as \code{p}.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing. This should be
#'  in the same units as the coordinates of \code{p} and \code{shoreline}.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing (see details).
#' @return A named vector representing the fetch length for each direction
#'  given in \code{bearings}.
#' @seealso \code{\link{fetch_len}} to work in geographic (unprojected)
#'  coordinates.
#' @export
fetch_len_proj <- function(p, bearings, shoreline, dmax, spread = 0) {
    # Check inputs
    if (!is(p, "SpatialPoints") || length(p) != 1) {
       stop("p must be a SpatialPoints* object of length 1.")
    }
    if (!(is(shoreline, "SpatialLines") || is(shoreline, "SpatialPolygons"))) {
        stop("shoreline must be a SpatialLines* or SpatialPolygons* object.")
    }
    if (proj4string(p) != proj4string(shoreline)) {
        stop("projections of p and shoreline do not match.")
    }
    if (!is.vector(bearings, "numeric")) stop("bearings must be a numeric vector.")
    if (!is.vector(spread, "numeric")) stop("spread must be a numeric vector.")
    if (!is.vector(dmax, "numeric") || length(dmax) != 1) {
        stop("dmax must be a single number.")
    }

    # If shoreline is a polygons (land) layer, check that point is not on land
    if (is(shoreline, "SpatialPolygons")) {
        # Note: 'as' only to remove DataFrame part of Spatial objects
        in_water <- is.na(over(as(p, "SpatialPoints"),
                               as(shoreline, "SpatialPolygons")))
        if(!in_water) {
            warning("point on land, returning NA")
            return(rep(NA, length(bearings)))
        }
    }

    # Clip shoreline layer to a rectangle around point
    # to guarantee at least dmax on each side
    clip_rect <- get_clip_rect_proj(p, dmax)
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
                            function(b) dist_shore_proj(p, shore_clip, b, dmax), 0)
    } else {
        # calculate the distance to shore for each sub-bearing
        bear_mat <- outer(bearings, spread, "+")
        dists <- vapply(bear_mat,
                        function(b) dist_shore_proj(p, shore_clip, b, dmax), 0)
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


# Helper functions below are not exported by the package

# Create clipping rectangle around point p
#   to guarantee at least dmax on each side
get_clip_rect_proj <- function(p, dmax) {
    x <- coordinates(p)[1]
    y <- coordinates(p)[2]
    r1 <- poly_rect(x - dmax, y - dmax, x + dmax, y + dmax)
    SpatialPolygons(list(Polygons(list(r1), ID = 1)),
                    proj4string = CRS(proj4string(p)))
}

# Returns the distance from point p to shoreline
# following bearing bear and up to distance dmax
dist_shore_proj <- function(p, shoreline, bear, dmax) {
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
}
