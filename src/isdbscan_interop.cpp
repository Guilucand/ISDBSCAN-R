
#include "structs.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

using namespace Rcpp;

IsdbscanResult isdbscan(std::vector<kd_point<double>> dataset, int k, int batch_size, bool stratif, bool approximate);

static Rcpp::List isdbscanToRcppList(IsdbscanResult const& result) {

    if (!result.layer.empty() && !result.border.empty()) {
        return Rcpp::List::create(
                Rcpp::Named("clusters") = result.clusters,
                Rcpp::Named("layer") = result.layer,
                Rcpp::Named("border") = result.border);
    }
    else {
        return Rcpp::List::create(
                Rcpp::Named("clusters") = result.clusters
        );
    }
}

template <class T>
std::vector<kd_point<double>> load_points(T const& data, bool transpose) {
    int nsamples = data->get_nrow();
    int nfeatures = data->get_ncol();

    if (transpose) std::swap(nsamples, nfeatures);

    std::vector<kd_point<double>> result;
    result.reserve(nsamples);

    Rcpp::NumericVector tmp_col(nfeatures);

    for (int i = 0; i < nsamples; i++) {
        if (transpose) {
            data->get_col(i, tmp_col.begin());
        } else {
            data->get_row(i, tmp_col.begin());
        }
        auto point = kd_point<double>();
        point.index = i;
        point.features.resize(nfeatures);
        for (int j = 0; j < nfeatures; j++) {
            point.features[j] = tmp_col(j);
        }
        result.push_back(point);
    }
    return result;
}

//'
//' isdbscan
//'
//'@param dataset numeric matrix or integer matrix or HDF5Matrix
//'@param k the number of clusters
//'@param batch_size the size of the mini batches
//'@details
//'This function performs the isdbscan algorithm over the specified dataset.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List isdbscan(SEXP data, int k, int batch_size, bool stratif, bool approximate, bool transpose) {
    auto data_rmatrix = beachmat::create_numeric_matrix(data);
    auto dataset = load_points(data_rmatrix, transpose);

    return isdbscanToRcppList(isdbscan(dataset, k, batch_size, stratif, approximate));

}