#include "rcpp_get_coocurrence_matrix.h"
#include "rcpp_create_neighborhood.h"
#include "rcpp_get_unique_values.h"
#include "get_class_index_map.h"

// [[Rcpp::export]]
IntegerMatrix rcpp_get_coocurrence_matrix(const IntegerMatrix x,
                                          const arma::imat directions,
                                          const int n_cores) {

    const int na = NA_INTEGER;
    const unsigned ncols = x.ncol();
    const unsigned nrows = x.nrow();

    std::vector<int> classes = rcpp_get_unique_values(x);
    std::map<int, unsigned> class_index = get_class_index_map(classes);

    const unsigned n_classes = class_index.size();
    std::vector<std::vector<unsigned> > cooc_mat(n_classes,
                                                 std::vector<unsigned>(n_classes));

    // create neighbors coordinates
    IntegerMatrix tmp = rcpp_create_neighborhood(directions);
    const int neigh_len = tmp.nrow();
    std::vector<std::vector<int> > neig_coords;
    for (int row = 0; row < neigh_len; row++) {
        IntegerVector a = tmp.row(row);
        const std::vector<int> b(a.begin(), a.end());
        neig_coords.push_back(b);
    }

    // NAs need an index, otherwise they are counted as neighbors of class[0]
    class_index.insert(std::make_pair(na, n_classes));

    // Setting the cores
    omp_set_num_threads(n_cores);
#pragma omp parallel default(none) shared(class_index, neig_coords, cooc_mat)
{
    // per thread
    std::vector<std::vector<unsigned> > cooc_mat_par(n_classes,
                                                     std::vector<unsigned>(n_classes));
#pragma omp for
    for (unsigned col = 0; col < ncols; col++) {
        for (unsigned row = 0; row < nrows; row++) {
            const int tmp = x[col * nrows + row];
            if (tmp == na)
                continue;
            unsigned focal_class = class_index[tmp];
            for (int h = 0; h < neigh_len; h++) {
                int neig_col = neig_coords[h][0] + col;
                int neig_row = neig_coords[h][1] + row;
                if (neig_col >= 0 &&
                    neig_row >= 0 &&
                    neig_col < ncols &&
                    neig_row < nrows) {
                    const int tmp = x[neig_col * nrows + neig_row];
                    if (tmp == na)
                        continue;
                    unsigned neig_class = class_index[tmp];
                    cooc_mat_par[focal_class][neig_class]++;
                }
            }
        }
    }

    // copy the parallel partial results to the final cooc_mat
    for(auto focal_class = 0; focal_class < n_classes; focal_class++) {
        for(auto neig_class = 0; neig_class < n_classes; neig_class++) {
#pragma omp critical
            cooc_mat[focal_class][neig_class] += cooc_mat_par[focal_class][neig_class];
        }
    }
}

IntegerMatrix result(n_classes, n_classes);
for (unsigned col = 0; col < cooc_mat.size(); col++) {
    for (unsigned row = 0; row < cooc_mat[col].size(); row++) {
        result(col, row) = cooc_mat[col][row];
    }
}

// add names
List u_names = List::create(classes, classes);
result.attr("dimnames") = u_names;
return result;
}

/*** R

library(raster)
library(dplyr)
test <- landscapemetrics::augusta_nlcd
mat <- raster::as.matrix(test)
four <- as.matrix(4)

test <- raster("~/Downloads/lc_2008_4bit_clip.tif")

not_par = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four)
cores_1 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 1)
cores_2 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 2)
cores_3 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 3)
cores_4 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 4)

identical(not_par, cores_4)

bench::mark(
    not_par = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four),
    cores_1 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 1),
    cores_2 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 2),
    cores_3 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 3),
    cores_4 = landscapemetrics:::rcpp_get_coocurrence_matrix(mat, four, 4),
    iterations = 1000,
    check = FALSE)

rcpp_get_coocurrence_matrix(mat, four)
rcpp_get_coocurrence_matrix(mat, four, 1)
rcpp_get_coocurrence_matrix(mat, four, 2)
rcpp_get_coocurrence_matrix(mat, four, 5)

test <- NLMR::nlm_mpd(5000, 5000) %>%
    landscapetools::util_classify(weighting = c(0.5, 0.15, 0.15, 0.05, 0.05, 0.1))

mat <- raster::as.matrix(test)

rcpp_get_coocurrence_matrix(mat, four)

Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

rcpp_get_coocurrence_matrix_par(mat, four, 5)

lsm_p_contig(test)

rcpp_get_unique_values(mat)
*/
