context("pad_raster")

test_that("pad_raster can handle all raster inputs", {

    landscape_matrix <- raster::as.matrix(landscape)

    expect_is(pad_raster(landscape), "list")
    expect_is(pad_raster(landscape_stack), "list")
    expect_is(pad_raster(landscape_brick), "list")
    expect_is(pad_raster(landscape_list), "list")
    expect_is(pad_raster(landscape_matrix, return_raster = FALSE, to_disk = FALSE), "list")
})

test_that("pad_raster adds correct number of rows/cols", {

    result_a <- landscapemetrics:::pad_raster(landscape,
                                              pad_raster_cells = 1,
                                              global = TRUE)[[1]]

    result_b <- landscapemetrics:::pad_raster(landscape,
                                              pad_raster_cells = 5,
                                              global = TRUE)[[1]]

    # original landscape is 30 x 30 and adding 1 rows/cols in each direction
    expect_equal(object = raster::nrow(result_a), expected = 32)
    expect_equal(object = raster::nrow(result_a), expected = 32)

    # original landscape is 30 x 30 and adding 5 rows/cols in each direction
    expect_equal(object = raster::nrow(result_b), expected = 40)
    expect_equal(object = raster::nrow(result_b), expected = 40)
})

test_that("pad_raster can add different padding values", {

    result <- landscapemetrics:::pad_raster(landscape,
                                            pad_raster_value = NA)[[1]]

    expect_true(anyNA(raster::values(result)))
})

test_that("pad_raster can return RasterLayer and matrix", {

    raster_layer <- landscapemetrics:::pad_raster(landscape)[[1]]

    mat <- landscapemetrics:::pad_raster(landscape, return_raster = FALSE)[[1]]

    expect_is(object = raster_layer, class = "RasterLayer")

    expect_is(object = mat, class = "matrix")
})
