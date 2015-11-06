library(waver)
context("Fetch length")

plonglat <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

p <- SpatialPoints(matrix(c(0, 1), ncol = 2), proj4string = CRS(plonglat))
q <- SpatialPoints(matrix(c(180, 25), ncol = 2), proj4string = CRS(plonglat))

bearings <- 90 * 0:3
dmax <- 1000
shoreline <- waver:::get_clip_rect(coordinates(p), dmax)


test_that("error when p not correctly specified for fetch_len", {
    expect_error(fetch_len(rbind(p, q), bearings, shoreline, dmax),
                 "single point")
    expect_error(fetch_len(0, bearings, shoreline, dmax), "single point")
    expect_error(fetch_len(cbind(c(1,2), c(2,0)), bearings, shoreline, dmax),
                 "single point")
})
