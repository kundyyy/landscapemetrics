#' RELMUTINF (landscape level)
#'
#' @description Relative mutual information
#'
#' @param landscape Raster* Layer, Stack, Brick or a list of rasterLayers.
#' @param neighbourhood The number of directions in which cell adjacencies are considered as neighbours:
#' 4 (rook's case) or 8 (queen's case). The default is 4.
#' @param ordered The type of pairs considered.
#' Either ordered (TRUE) or unordered (FALSE).
#' The default is TRUE.
#' @param base The unit in which entropy is measured.
#' The default is "log2", which compute entropy in "bits".
#' "log" and "log10" can be also used.
#'
#' @details
#' Due to the spatial autocorrelation, the value of mutual information tends to grow with a diversity of the landscape (marginal entropy). To adjust this tendency, it is possible to calculate relative mutual information by dividing the mutual information by the marginal entropy. Relative mutual information always has a range between 0 and 1 and can be used to compare spatial data with different number and distribution of categories.
#' When the value of mutual information equals to 0, then relative mutual information is 1.
#'
#' @seealso
#' \code{\link{lsm_l_ent}},
#' \code{\link{lsm_l_condent}},
#' \code{\link{lsm_l_joinent}},
#' \code{\link{lsm_l_mutinf}}
#'
#' @return tibble
#'
#' @examples
#' lsm_l_relmutinf(landscape)
#'
#' @aliases lsm_l_relmutinf
#' @rdname lsm_l_relmutinf
#'
#' @references
#' Nowosad J., TF Stepinski. 2019. Information theory as a consistent framework
#' for quantification and classification of landscape patterns. https://doi.org/10.1007/s10980-019-00830-x
#'
#' @export
lsm_l_relmutinf <- function(landscape,
                         neighbourhood,
                         ordered,
                         base) UseMethod("lsm_l_relmutinf")

#' @name lsm_l_relmutinf
#' @export
lsm_l_relmutinf.RasterLayer <- function(landscape,
                                     neighbourhood = 4,
                                     ordered = TRUE,
                                     base = "log2") {

    result <- lapply(X = raster::as.list(landscape),
                     FUN = lsm_l_relmutinf_calc,
                     neighbourhood = neighbourhood,
                     ordered = ordered,
                     base = base)

    layer <- rep(seq_along(result),
                 vapply(result, nrow, FUN.VALUE = integer(1)))

    result <- do.call(rbind, result)

    tibble::add_column(result, layer, .before = TRUE)
}

#' @name lsm_l_relmutinf
#' @export
lsm_l_relmutinf.RasterStack <- function(landscape,
                                     neighbourhood = 4,
                                     ordered = TRUE,
                                     base = "log2") {

    result <- lapply(X = raster::as.list(landscape),
                     FUN = lsm_l_relmutinf_calc,
                     neighbourhood = neighbourhood,
                     ordered = ordered,
                     base = base)

    layer <- rep(seq_along(result),
                 vapply(result, nrow, FUN.VALUE = integer(1)))

    result <- do.call(rbind, result)

    tibble::add_column(result, layer, .before = TRUE)
}

#' @name lsm_l_relmutinf
#' @export
lsm_l_relmutinf.RasterBrick <- function(landscape,
                                     neighbourhood = 4,
                                     ordered = TRUE,
                                     base = "log2") {

    result <- lapply(X = raster::as.list(landscape),
                     FUN = lsm_l_relmutinf_calc,
                     neighbourhood = neighbourhood,
                     ordered = ordered,
                     base = base)

    layer <- rep(seq_along(result),
                 vapply(result, nrow, FUN.VALUE = integer(1)))

    result <- do.call(rbind, result)

    tibble::add_column(result, layer, .before = TRUE)
}

#' @name lsm_l_relmutinf
#' @export
lsm_l_relmutinf.stars <- function(landscape,
                                neighbourhood = 4,
                                ordered = TRUE,
                                base = "log2") {

    landscape <- methods::as(landscape, "Raster")

    result <- lapply(X = raster::as.list(landscape),
                     FUN = lsm_l_relmutinf_calc,
                     neighbourhood = neighbourhood,
                     ordered = ordered,
                     base = base)

    layer <- rep(seq_along(result),
                 vapply(result, nrow, FUN.VALUE = integer(1)))

    result <- do.call(rbind, result)

    tibble::add_column(result, layer, .before = TRUE)
}

#' @name lsm_l_relmutinf
#' @export
lsm_l_relmutinf.list <- function(landscape,
                              neighbourhood = 4,
                              ordered = TRUE,
                              base = "log2") {

    result <- lapply(X = landscape,
                     FUN = lsm_l_relmutinf_calc,
                     neighbourhood = neighbourhood,
                     ordered = ordered,
                     base = base)

    layer <- rep(seq_along(result),
                 vapply(result, nrow, FUN.VALUE = integer(1)))

    result <- do.call(rbind, result)

    tibble::add_column(result, layer, .before = TRUE)
}

lsm_l_relmutinf_calc <- function(landscape, neighbourhood, ordered, base){

    # convert to matrix
    if (!inherits(x = landscape, what = "matrix")) {
        landscape <- raster::as.matrix(landscape)
    }

    # all values NA
    if (all(is.na(landscape))) {
        return(tibble::tibble(level = "landscape",
                              class = as.integer(NA),
                              id = as.integer(NA),
                              metric = "mutinf",
                              value = as.double(NA)))
    }

    com <- rcpp_get_coocurrence_matrix(landscape,
                                       directions = as.matrix(neighbourhood))
    com_c <- colSums(com)

    coh <- rcpp_get_coocurrence_vector(landscape,
                                       directions = as.matrix(neighbourhood),
                                       ordered = ordered)

    comp <- rcpp_get_entropy(com_c, base)
    cplx <- rcpp_get_entropy(coh, base)
    conf <- cplx - comp
    aggr <- comp - conf
    rel <- ifelse(aggr == 0, 1, aggr / comp)

    return(tibble::tibble(level = "landscape",
                          class = as.integer(NA),
                          id = as.integer(NA),
                          metric = "relmutinf",
                          value = as.double(rel)))
}
