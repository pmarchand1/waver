library(waver)
context("Fetch length")

sf_use_s2(FALSE)

# Set up Spatial objects
longlat <- st_crs("EPSG:4326")
aeqd1 <- st_crs("ESRI:54032")
aeqd2 <- st_crs("+proj=aeqd +lon_0=180 +lat_0=0")

# p1 at (0,0)
p1 <- st_sfc(st_point(c(0, 0)), crs = longlat)
land1 <- st_sfc(waver:::poly_rect(-0.2, 0.25, 0.3, 0.5), crs = longlat)
p1_df <- st_sf(v1 = 2, v2 = "a", p1)
p1_prj <- st_transform(p1, aeqd1)
land1_prj <- st_transform(land1, aeqd1)

# p2 near international date line
p2 <- st_sfc(st_point(c(-179.75, 0)), crs = longlat)
land2 <- st_sfc(list(
    waver:::poly_rect(-179.95, 0.25, -179.45, 0.5),
    waver:::poly_rect(179.7, -0.4, 179.95, 0.2)
), crs = longlat)
p2_prj <- st_transform(p2, aeqd2)
land2_prj <- st_transform(land2, aeqd2)

# p3 on land
p3 <- st_sfc(st_point(c(0, 0.4)), crs = longlat)
p3_prj <- st_transform(p3, aeqd1)
p13 <- c(p1, p3)
p13_prj <- st_transform(p13, aeqd1)

# Multiple points
pts <- c(p1, p2, p3)
lands <- c(land1, land2)
pts_df <- st_sf(v1 = c(2, 5, 6), v2 = c("a", "e", "h"), pts)
lands_df <- st_sf(v3 = c("s", "t", "u"), lands)

# Fetch parameters
bearings <- c(0, 45, 225, 315)
spread <- c(-10, 0, 10)
dmax <- 50000

# Expected results
tol <- 1E-4 # relative tolerance level
fexp1 <- setNames(c(27644, 39094, 50000, 50000), bearings)
fexp1spr <- setNames(c(27926, 40937, 50000, 44610), bearings)
fexp2 <- setNames(c(27644, 39094, 47228, 50000), bearings)
fexp2spr <- setNames(c(27926, 40937, 46005, 44610), bearings)
fexp3 <- setNames(rep(NA, 4), bearings)


test_that("fetch_len correct for single point", {
    expect_equal(fetch_len(p1, bearings, land1, dmax), fexp1, tolerance = tol)
    expect_equal(fetch_len(p1, bearings, land1, dmax, spread),
                 fexp1spr, tolerance = tol)
    expect_equal(fetch_len(p1, bearings, st_cast(land1, "LINESTRING"), dmax),
                 fexp1, tolerance = tol)
    expect_equal(fetch_len(p1_prj, bearings, land1_prj, dmax),
                 fexp1, tolerance = tol)
    expect_equal(fetch_len(p1_prj, bearings, land1_prj, dmax, spread),
                 fexp1spr, tolerance = tol)
    expect_equal(fetch_len(p1_prj, bearings, st_cast(land1_prj, "LINESTRING"), dmax),
                 fexp1, tolerance = tol)
})


test_that("fetch_len correct near international date line", {
    expect_equal(fetch_len(p2, bearings, land2, dmax), fexp2, tolerance = tol)
    expect_equal(fetch_len(p2, bearings, land2, dmax, spread),
                 fexp2spr, tolerance = tol)
    expect_equal(fetch_len(p2_prj, bearings, land2_prj, dmax),
                 fexp2, tolerance = tol)
    expect_equal(fetch_len(p2_prj, bearings, land2_prj, dmax, spread),
                 fexp2spr, tolerance = tol)
})


test_that("fetch_len for point on land returns NA and issues warning", {
    expect_warning(fetch_len(p3, bearings, land1, dmax), "on land")
    expect_warning(fetch_len(p3_prj, bearings, land1_prj, dmax), "on land")
    expect_equal(suppressWarnings(fetch_len(p3, bearings, land1, dmax)), fexp3)
    expect_equal(suppressWarnings(fetch_len(p3_prj, bearings, land1_prj, dmax)), fexp3)
})


test_that("fetch_len_multi matches fetch_len results", {
    expect_equal(suppressWarnings(fetch_len_multi(pts, bearings, lands, dmax)),
                 `rownames<-`(rbind(fexp1, fexp2, fexp3), NULL), tolerance = tol)
    expect_equal(suppressWarnings(fetch_len_multi(pts, bearings, lands, dmax,
                                                  method = "clip")),
                 `rownames<-`(rbind(fexp1, fexp2, fexp3), NULL), tolerance = tol)
    expect_equal(suppressWarnings(fetch_len_multi(pts, bearings, lands, dmax,
                                                  method = "clip")),
                 `rownames<-`(rbind(fexp1, fexp2, fexp3), NULL), tolerance = tol)
    expect_equal(
        suppressWarnings(fetch_len_multi(pts, bearings, lands, dmax, spread)),
        `rownames<-`(rbind(fexp1spr, fexp2spr, fexp3), NULL), tolerance = tol)
    expect_equal(
        suppressWarnings(fetch_len_multi(p13_prj, bearings, land1_prj, dmax)),
        `rownames<-`(rbind(fexp1, fexp3), NULL), tolerance = tol)
})


test_that("fetch_len works with sf objects", {
    expect_equal(fetch_len(p1_df, bearings, lands_df, dmax),
                 fexp1, tolerance = tol)
    expect_equal(
        suppressWarnings(fetch_len_multi(pts_df, bearings, lands_df, dmax)),
        `rownames<-`(rbind(fexp1, fexp2, fexp3), NULL), tolerance = tol)
})


test_that("fetch_len fails on bad inputs", {
    expect_error(fetch_len(0, bearings, land1, dmax), "spatial point")
    expect_error(fetch_len(pts, bearings, land1, dmax), "single point")
    expect_error(fetch_len(p1, bearings, p1, dmax), "lines or polygons")
    expect_error(fetch_len(p1_prj, bearings, land1, dmax), "projections")
    expect_error(fetch_len(p1, bearings, land1_prj, dmax), "projections")
    expect_error(fetch_len(p1_prj, bearings, st_transform(land1, aeqd2),
                           dmax), "projections")
    expect_error(fetch_len(p1, "a", land1, dmax), "bearings")
    expect_error(fetch_len(p1, bearings, land1, dmax, "a"), "spread")
    expect_error(fetch_len(p1, bearings, land1, 0), "dmax")
    expect_error(fetch_len(p1, bearings, land1, c(1,2)), "dmax")
})
