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
bearing_line <- function(p, bearing, len) {
    l1 <- Line(rbind(
        coordinates(p),
        coordinates(p) + len * c(sinpi(bearing / 180), cospi(bearing / 180))
    ))
    SpatialLines(list(Lines(list(l1), ID = 1)), CRS(proj4string(p)))
}
