#### Some utility functions to create / manipulate spatial objects (not exported)

# Create a polygon (simple feature) corresponding to rectangle with given coords
poly_rect <- function(xmin, ymin, xmax, ymax) {
    st_polygon(list(cbind(c(rep(xmin, 2), rep(xmax, 2), xmin),
                          c(ymin, rep(ymax, 2), rep(ymin, 2)))))
}

# Create a linestring (simple feature) starting at p, of length len,
#     along given bearing (in degrees)
# Note: This function is meant to be used with a projected CRS
bearing_line <- function(p, bearing, len) {
    l1 <- st_linestring(rbind(
        st_coordinates(p),
        st_coordinates(p) + len * c(sinpi(bearing / 180), cospi(bearing / 180))
    ))
    st_sfc(l1, crs = st_crs(p))
}

# Create clipping rectangle around point p (sfc_POINT object for a single point)
#  to guarantee at least dmax on each side
# dmax either in meters (if longlat coordinates) or in the projection's coordinates
get_clip_rect <- function(p, dmax) {
    if (st_is_longlat(p)) {
        lat_dist <- 111600 # approx. distance (in m) between degrees of latitude
        long <- st_coordinates(p)[1]
        lat <- st_coordinates(p)[2]
        ybuf <- dmax / lat_dist
        xbuf <- ybuf / cospi(abs(lat) / 180)
        # Split clip_rect in two if it would overlap international date line
        if (long - xbuf < -180) {
            westr <- poly_rect(-180, lat - ybuf, long + xbuf, lat + ybuf)
            eastr <- poly_rect(long - xbuf + 360, lat - ybuf, 180, lat + ybuf)
            suppressMessages(
                clip_rect <- st_union(st_sfc(westr, eastr, crs = st_crs(p)))
            )
        } else if(long + xbuf > 180) {
            westr <- poly_rect(-180, lat - ybuf, long + xbuf - 360, lat + ybuf)
            eastr <- poly_rect(long - xbuf, lat - ybuf, 180, lat + ybuf)
            suppressMessages(
                clip_rect <- st_union(st_sfc(westr, eastr, crs = st_crs(p)))
            )
        } else {
            r1 <- poly_rect(long - xbuf, lat - ybuf, long + xbuf, lat + ybuf)
            clip_rect <- st_sfc(r1, crs = st_crs(p))
        }
    } else {
        x <- st_coordinates(p)[1]
        y <- st_coordinates(p)[2]
        r1 <- poly_rect(x - dmax, y - dmax, x + dmax, y + dmax)
        clip_rect <- st_sfc(r1, crs = st_crs(p))
    }
    clip_rect
}

