# Some utility functions to create / manipulate spatial objects


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
