// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

#include <iostream>
#include <algorithm>
#include "structs.h"
#include "knn.h"
#include "influence_space.h"
#include "clustering.h"
#include "inflo.h"
#include "stratifier.h"

static void print_frequencies(std::vector<knn_influence_space> const& influence_space, int k) {

    std::vector<int> cnts(k + 1, 0);
    for (auto &x : influence_space) {
        cnts[x.size()]++;
    }

    std::cout << "influence_space FREQUENCIES: " << std::endl;
    for (int i = 0; i < k + 1; i++) {
        std::cout << "ksize = " << i << " count = " << cnts[i] << std::endl;
    }
}

std::vector<double> calculate_nnk_dist(std::vector<knn_dists<double>> const& knn_dists, int k, double &max_nndist){

    std::vector<double> nnk_dist(knn_dists.size(), 0);

    for(int i = 0; i < knn_dists.size(); i++){
        double sum = 0;
        for(int j = 0; j < k; j++){
            sum += knn_dists[i][j].distance;
        }
        nnk_dist[i] = sum;
        max_nndist = std::max(max_nndist, sum);
    }
}

std::vector<int> calculate_rnn(std::vector<knn_dists<double>> const& knn_dists, int k){

    std::vector<int> nof_rnn(knn_dists.size());

    for(const auto& knn_dist : knn_dists) {
        for (int j = 0; j < k; j++) {
            nof_rnn[knn_dist[j].index]++;
        }
    }
}

std::vector<std::pair<double, int>> sort_points(std::vector<double> const& inflo, std::vector<double> const& nnk_dists) {

    std::vector<std::pair<double, int>> result(inflo.size());

    for (int i = 0; i < result.size(); i++) {
        result[i] = std::make_pair(i, inflo[i] + nnk_dists[i]);
    }

    std::sort(result.begin(), result.end());

    return result;
}

std::vector<int> build_ranks_map(std::vector<std::pair<double, int>> const& sorted_pts) {
    std::vector<int> result(sorted_pts.size());
    for (int i = 0; i < result.size(); i++) {
        result[sorted_pts[i].second] = i;
    }
    return result;
}

IsdbscanResult isdbscan(std::vector<kd_point<double>> dataset, int k, int batch_size, bool stratif, bool approximate) {

    std::cout << "Dataset size: (" << dataset.size() << "x" << (dataset.size() ? dataset.front().features.size() : 0) << ")" << std::endl;

    std::vector<int> rnn_count(dataset.size());
    std::vector<bool> border(dataset.size());

    // K cannot be larger than the dataset - 1
    if (k > (int)dataset.size() - 1) {
        k = (int) dataset.size() - 1;
        std::cout << "WARNING: K set to " << k << std::endl;
    }

    //When data fits memory
    if (dataset.size() <= batch_size) {

        auto knn_dists = calc_knn_dists(dataset, k, approximate);

        auto influence_space = find_influence_space(knn_dists, k);
        print_frequencies(influence_space, k);

        if (!stratif){
            auto clustering_res = calc_cluster(influence_space, k, border);
            return IsdbscanResult {
                /*.clusters = */ clustering_res
            };
        }
//    Rcpp::NumericMatrix submat(data_n_rows, data_n_cols);
//    Rcpp::IntegerMatrix calIsk(data_n_rows,data_n_cols);
//    Rcpp::LogicalVector border(data_n_rows);
        else{
            double max_inflo;
            double max_nndist;

            auto nnk_dist = calculate_nnk_dist(knn_dists, k, max_nndist);
            auto nof_rnn = calculate_rnn(knn_dists, k);
            std::vector<double> inflo_vector =
                    calculate_inflo(knn_dists, influence_space, max_inflo);
            calculate_inflo_correct(nof_rnn, influence_space, inflo_vector);
            calculate_inflo_norm(inflo_vector, max_inflo, max_nndist);

            auto inflo_plus_nn_ordered = sort_points(inflo_vector, nnk_dist);
            auto ranks_map = build_ranks_map(inflo_plus_nn_ordered);

            auto layer_vec = stratifier(inflo_vector, inflo_plus_nn_ordered, ranks_map);
//
            auto clustering_res = calc_cluster(influence_space, k, border);
            return IsdbscanResult {
                    /*.clusters = */ clustering_res,
                    /* .layer = */ layer_vec,
                    /* .border = */ border
            };
        }
    }
//    else{
//
//        auto final_matrix=beachmat::create_numeric_matrix(data);
//        int rowsToRead = std::floor(batch_size/2);
//        int howManyTimesTofitData = std::ceil(data_n_rows/rowsToRead);
//        Rcpp::NumericVector tmpForj(data_n_cols);
//        int offset;
//        int counter = 0;
//
//        Rcpp::IntegerVector pointMaxIdx(data_n_rows);
//        Rcpp::NumericVector pointMaxValue(data_n_rows);
//
//        for (int i = 0; i < data_n_rows; i++){
//            pointMaxIdx(i) = -1;
//            pointMaxValue(i) = 99999;
//        }
//
//        for (int i = 0; i < data_n_rows; i++){
//            final_matrix->get_row(i, tmp.begin());
//            for (counter = 0; counter < howManyTimesTofitData; counter++){
//                offset = rowsToRead * counter;
//                for(int j = rowsToRead * counter;  j < data_n_rows; j++){
//                    final_matrix->get_row(j, tmpForj.begin());
//                    submat.row(j-offset) = tmpForj;
//                }
//                calcKnnDists(tmp ,submat, knnDists, knnIdx, pointMaxIdx, pointMaxValue, k, data_n_rows, data_n_cols, offset, rowsToRead, i);
//            }
//        }
//
//        sortKnnDists(knnDists, knnIdx, k, data_n_rows);
//        calIsk = findIsk(knnIdx,k,data_n_rows,nOfIsk);
//        if (!stratif){
//            Rcpp::IntegerVector clusteringRes = calc_cluster(calIsk,data_n_rows,k,nOfIsk,border);
//            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes);
//        }
//        else{
//            double max_NNK = -1;
//            double max_inflo = -1;
//            double max_NNdist = -1;
//            Rcpp::NumericVector NNdistVector(data_n_rows); //da creare solo se stratification == true;
//            Rcpp::NumericVector INFLOVector(data_n_rows);
//            Rcpp::IntegerVector map(data_n_rows);
//            Rcpp::IntegerVector newOrder(data_n_rows);
//            Rcpp::NumericVector INFLOPlusNNdist(data_n_rows);
//            Rcpp::IntegerVector layerVector(data_n_rows);
//
//
//            calculateNNKdist(knnDists,NNdistVector,k,data_n_rows, max_NNdist);
//            calculateRNN(knnIdx,nOfRNN,k,data_n_rows);
//            calculateINFLO(knnDists,knnIdx,nOfIsk,INFLOVector, k, data_n_rows, max_inflo);
//            calculate_inflo_correct(nOfRNN, nOfIsk, INFLOVector, data_n_rows);
//            calculateINFLONorm(INFLOVector, data_n_rows, max_inflo, max_NNdist);
//            sortPoints(INFLOVector,NNdistVector, map,newOrder, INFLOPlusNNdist, data_n_rows);
//            stratifier(INFLOVector,NNdistVector,INFLOPlusNNdist,map,newOrder,layerVector,data_n_rows);
//
//            Rcpp::IntegerVector clusteringRes = calc_cluster(calIsk,data_n_rows,k,nOfIsk,border,map);
//            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes, Rcpp::Named("layer") = layerVector, Rcpp::Named("border") = border);
//        }
//    }
}