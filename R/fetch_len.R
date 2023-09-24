
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
#' The input data can be in either geographic (long, lat) or projected coordinates,
#' but \code{p} and \code{shoreline} must share the same coordinate system. Distances
#' are calculated using the \code{\link[sf]{st_distance}} function from the sf package
#' and expressed in the units of the coordinate system used, or in meters if using
#' geographic coordinates. For geographic coordinates, we recommend setting
#' \code{sf_use_s2(FALSE)}, which results in \code{st_distance} using the ellipsoid
#' distance calculation (requires the lwgeom package), instead of the less precise
#'  spherical distance calculation. For projected coordinates, the Euclidean distance
#' is calculated.
#'
#' If the shoreline layer is composed of polygons rather than lines, the function
#' verifies that the input point is outside all polygons (i.e. in water). If this is
#' not the case, it issues a warning and returns a vector of \code{NA}.
#'
#' @param p Simple feature (sf or sfc) object representing a single point.
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline Simple feature (sf or sfc) object representing the
#'  shoreline, in either line or polygon format.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing (see details).
#' @param projected Deprecated argument, kept for backwards compatibility.
#' @param check_inputs Should the validity of inputs be checked? It is
#' recommended to keep this TRUE, unless this function is called repeatedly from
#' another function that already checks inputs.
#' @return A named vector representing the fetch length for each direction
#'  given in \code{bearings}.
#' @examples
#'  pt <- st_sfc(st_point(c(0, 0)), crs = st_crs(4326))
#'  # Shoreline is a rectangle from (-0.2, 0.25) to (0.3, 0.5)
#'  rect <- st_polygon(list(cbind(c(rep(-0.2, 2), rep(0.3, 2), -0.2),
#'                                c(0.25, rep(0.3, 2), rep(0.25, 2)))))
#'  land <- st_sfc(rect, crs = st_crs(4326))
#'  fetch_len(pt, bearings = c(0, 45, 225, 315), land,
#'            dmax = 50000, spread = c(-10, 0, 10))
#' @seealso \code{\link{fetch_len_multi}} for an efficient alternative when
#'  computing fetch length for multiple points.
#' @export
fetch_len <- function(p, bearings, shoreline, dmax, spread = 0,
                      projected = FALSE, check_inputs = TRUE) {
    # Convert sp objects to sf format
    if(is(p, "SpatialPointsDataFrame")) p <- st_as_sf(p)
    if(is(p, "SpatialPoints")) p <- st_as_sfc(p)
    if(is(shoreline, "SpatialLinesDataFrame") || is(shoreline, "SpatialPolygonsDataFrame")) {
        shoreline <- st_as_sf(shoreline)
    }
    if(is(shoreline, "SpatialLines") || is(shoreline, "SpatialPolygons")) {
        shoreline <- st_as_sfc(shoreline)
    }
    # Extract geometry columns if p or shoreline are sf objects
    if(is(p, "sf")) p <- st_geometry(p)
    if(is(shoreline, "sf")) shoreline <- st_geometry(shoreline)

    if (check_inputs) {
        if (!is(p, "sfc_POINT")) stop("p must be a spatial point in sf or sp format.")
        if(length(p) != 1) stop("p must be a single point.")
        if (!(is(shoreline, "sfc_LINESTRING") || is(shoreline, "sfc_MULTILINESTRING") ||
              is(shoreline, "sfc_POLYGON")    || is(shoreline, "sfc_MULTIPOLYGON"))) {
            stop("shoreline must be a spatial object in sf or sp format containing lines or polygons.")
        }
        if (st_crs(p) != st_crs(shoreline)) {
            stop("projections of p and shoreline do not match.")
        }

        if (!is.vector(bearings, "numeric")) stop("bearings must be a numeric vector.")
        if (!is.vector(spread, "numeric")) stop("spread must be a numeric vector.")
        if (!is.vector(dmax, "numeric") || length(dmax) != 1 || dmax <= 0) {
            stop("dmax must be a single number greater than 0.")
        }
    }

    # If shoreline is a polygons (land) layer, check that point is not on land
    if (is(shoreline, "sfc_POLYGON") || is(shoreline, "sfc_MULTIPOLYGON")) {
        suppressMessages(
            in_water <- length(st_intersection(p, shoreline)) == 0
        )
        if(!in_water) {
            warning("point on land, returning NA")
            return(setNames(rep(NA, length(bearings)), bearings))
        }
    }

    # Clip shoreline layer to a rectangle around point
    # to guarantee at least dmax on each side
    clip_rect <- get_clip_rect(p, dmax)
    suppressMessages(
        shore_clip <- st_intersection(shoreline, clip_rect)
    )

    # If no land within rectangle, return dmax for all bearings
    if (length(shore_clip) == 0) {
        return(setNames(rep(dmax, length(bearings)), bearings))
    }

    # Convert to multilinestring (if not already lines)
    #  since line-line intersections are needed later
    shore_clip <- st_cast(shore_clip, "MULTILINESTRING")

    # Calculate fetch
    if (all(spread == 0)) {
        # If no sub-bearings, just return distance to shore for each bearing
        fetch_res <- vapply(bearings,
                   function(b) dist_shore(p, shore_clip, b, dmax), 0)
    } else {
        # Calculate the distance to shore for each sub-bearing
        bear_mat <- outer(bearings, spread, "+")
        dists <- vapply(bear_mat,
                   function(b) dist_shore(p, shore_clip, b, dmax), 0)
        dim(dists) <- dim(bear_mat)
        # Return weighted means of the sub-bearing fetch values
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
#' With \code{method = "btree"} (default), the fetch calculation for each point only uses
#' the geometries within the \code{shoreline} layer that intersect with a rectangular
#' buffer of size \code{dmax} around that point. (The name is based on a previous version
#' of the function that implemented this method using the \code{gBinarySTRtreeQuery} function
#' from the rgeos package.)
#'
#' With \code{method = "clip"}, the \code{shoreline} is clipped to its intersection
#' with a polygon formed by the union of all the individual points' rectangular buffers.
#'
#' In both cases, \code{\link{fetch_len}} is then applied to each point,
#' using only the necessary portion of the shoreline.
#'
#' Generally, the "clip" method will produce the biggest time savings when
#' points are clustered within distances less than \code{dmax} (so their
#' clipping rectangles overlap), whereas the "btree" method will be more
#' efficient when the shoreline is composed of multiple geometrical objects
#' and points are distant from each other.
#'
#' @param pts Simple features (sf or sfc) object containing point data.
#' @param bearings Vector of bearings, in degrees.
#' @param shoreline Simple feature (sf or sfc) object representing the
#'  shoreline, in either line or polygon format.
#' @param dmax Maximum value of fetch length, returned if there is no land
#'  within a distance of \code{dmax} from a given bearing.
#' @param spread Vector of relative bearings (in degrees) for which
#'  to calculate fetch around each main bearing.
#' @param method Whether to use the "btree" (default) or "clip" method.
#'  See below for more details.
#' @param projected Deprecated argument, kept for backwards compatibility.
#' @return A matrix of fetch lengths, with one row by point in \code{pts} and
#'  one column by bearing in \code{bearings}.
#' @seealso \code{\link{fetch_len}} for details on the fetch length computation.
#' @export
fetch_len_multi <- function(pts, bearings, shoreline, dmax,
                     spread = 0, method = "btree", projected = FALSE) {
    # Convert sp objects to sf format
    if(is(pts, "SpatialPointsDataFrame")) pts <- st_as_sf(pts)
    if(is(pts, "SpatialPoints")) pts <- st_as_sfc(pts)
    if(is(shoreline, "SpatialLinesDataFrame") || is(shoreline, "SpatialPolygonsDataFrame")) {
        shoreline <- st_as_sf(shoreline)
    }
    if(is(shoreline, "SpatialLines") || is(shoreline, "SpatialPolygons")) {
        shoreline <- st_as_sfc(shoreline)
    }
    # Extract geometry columns if pts or shoreline are sf objects
    if(is(pts, "sf")) pts <- st_geometry(pts)
    if(is(shoreline, "sf")) shoreline <- st_geometry(shoreline)

    # Check inputs
    match.arg(method, choices = c("btree", "clip"))
    if (!is(pts, "sfc_POINT")) stop("pts must be a spatial point in sf or sp format.")
    if (!(is(shoreline, "sfc_LINESTRING") || is(shoreline, "sfc_MULTILINESTRING") ||
          is(shoreline, "sfc_POLYGON")    || is(shoreline, "sfc_MULTIPOLYGON"))) {
        stop("shoreline must be a spatial object in sf or sp format containing lines or polygons.")
    }
    if (st_crs(pts) != st_crs(shoreline)) {
        stop("projections of p and shoreline do not match.")
    }

    if (!is.vector(bearings, "numeric")) stop("bearings must be a numeric vector.")
    if (!is.vector(spread, "numeric")) stop("spread must be a numeric vector.")
    if (!is.vector(dmax, "numeric") || length(dmax) != 1 || dmax <= 0) {
        stop("dmax must be a single number greater than 0.")
    }

    # Create rectangular buffers around each point
    rect_list <- lapply(1:length(pts),
                        function(i) get_clip_rect(pts[i], dmax))
    rect_buf <- do.call(c, rect_list)

    if (method == "btree") {
        # Generate list of shoreline geometry indices that intersect each rectangle
        suppressMessages(
            btree <- st_intersects(rect_buf, shoreline)
        )
        # Calculate fetch for point at index i using btree
        fetch_i <- function(i) {
            if (length(btree[[i]]) == 0) {
                setNames(rep(dmax, length(bearings)), bearings)
            } else {
                fetch_len(pts[i], bearings, shoreline[btree[[i]]], dmax,
                          spread, check_inputs = FALSE)
            }
        }
        # Calculate fetch for all points and return a (points x bearings) matrix
        fetch_res <- t(vapply(1:length(pts), fetch_i, rep(0, length(bearings))))
    } else { # method == "clip"
        # Clip shoreline to a merged buffer around all points
        suppressMessages({
            rect_buf <- st_union(rect_buf)
            sub_shore <- st_intersection(shoreline, rect_buf)
        })
        fetch_res <- t(
            vapply(1:length(pts),
                   function(i) fetch_len(pts[i], bearings, sub_shore, dmax,
                                         spread, check_inputs = FALSE),
                   rep(0, length(bearings)))
        )
    }
    fetch_res
}


#### Helper functions below are not exported by the package ####

# Returns the distance from point p to shoreline following bearing bear and up to distance dmax
# Uses the geodetic distance if in geographic coordinates,
#   or the Euclidean distance if in projected coordinates
dist_shore <- function(p, shoreline, bear, dmax) {
    if (st_is_longlat(p)) {
        # Draw geodesic line of length dmax with given start point and bearing
        # Line drawn with a point every dmax/500
        geo_line <- geosphere::gcIntermediate(st_coordinates(p),
                                              geosphere::destPoint(st_coordinates(p), bear, dmax),
                                              n = 500, breakAtDateLine = TRUE, addStartEnd = TRUE
        )
        if (is.list(geo_line)) {
            geo_line <- st_multilinestring(geo_line)
        } else {
            geo_line <- st_linestring(geo_line)
        }
        bline <- st_sfc(geo_line, crs = st_crs(shoreline))
    } else {
        # Draw line of length dmax with given start point and bearing
        bline <- bearing_line(p, bear, dmax)
    }
    # Return (minimum) distance from p to intersection of bearing line and shoreline
    # If no intersection, fetch is dmax
    suppressMessages(
        land_int <- st_intersection(bline, shoreline)
    )
    if (length(land_int) == 0) {
        d <- dmax
    } else {
        d <- min(st_distance(p, land_int))
    }
    d
}
