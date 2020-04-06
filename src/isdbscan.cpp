// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

#include <iostream>
#include "structs.h"
#include "knn.h"
#include "biconnected.h"
#include "clustering.h"

static void print_frequencies(std::vector<knn_biconnected> const& biconnected, int k) {

    std::vector<int> cnts(k + 1, 0);
    for (auto &x : biconnected) {
        cnts[x.size()]++;
    }

    std::cout << "BICONNECTED FREQUENCIES: " << std::endl;
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

        auto biconnected = find_biconnected(knn_dists, k);
        print_frequencies(biconnected, k);

        if (!stratif){
            auto clusteringRes = calc_cluster(biconnected, k, border);
            return IsdbscanResult {
                /*.clusters = */ clusteringRes
            };
        }

        //    Rcpp::IntegerVector nOfIsk(data_n_rows);
//    Rcpp::NumericMatrix submat(data_n_rows, data_n_cols);
//    Rcpp::NumericMatrix knnDists(data_n_rows, k);
//    Rcpp::IntegerMatrix knnIdx(data_n_rows, k);
//    Rcpp::NumericVector tmp(data_n_cols);
//    Rcpp::IntegerVector nOfRNN(data_n_rows);
//    Rcpp::IntegerMatrix calIsk(data_n_rows,data_n_cols);
//    Rcpp::LogicalVector border(data_n_rows);
        else{
            double max_inflo;
            double max_nndist;
//            Rcpp::NumericVector NNdistVector(data_n_rows); //da creare solo se stratification == true;
//            Rcpp::NumericVector INFLOVector(data_n_rows);
//            Rcpp::IntegerVector map(data_n_rows);
//            Rcpp::IntegerVector newOrder(data_n_rows);
//            Rcpp::NumericVector INFLOPlusNNdist(data_n_rows);
//            Rcpp::IntegerVector layerVector(data_n_rows);
            auto nndist = calculate_nnk_dist(knn_dists, k, max_nndist);
            auto nof_rnn = calculate_rnn(knn_dists, k);
//            calculateINFLO(knnDists,knnIdx,nOfIsk,INFLOVector, k, data_n_rows, max_inflo);
//            calculateINFLOCorrect(nOfRNN, nOfIsk, INFLOVector, data_n_rows);
//            calculateINFLONorm(INFLOVector, data_n_rows, max_inflo, max_NNdist);
//            sortPoints(INFLOVector,NNdistVector, map,newOrder, INFLOPlusNNdist, data_n_rows);
//            stratifier(INFLOVector,NNdistVector,INFLOPlusNNdist,map,newOrder,layerVector,data_n_rows);
//
//            Rcpp::IntegerVector clusteringRes = calc_cluster(calIsk,data_n_rows,k,nOfIsk,border,map);
//            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes, Rcpp::Named("layer") = layerVector, Rcpp::Named("border") = border);
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
//            calculateINFLOCorrect(nOfRNN, nOfIsk, INFLOVector, data_n_rows);
//            calculateINFLONorm(INFLOVector, data_n_rows, max_inflo, max_NNdist);
//            sortPoints(INFLOVector,NNdistVector, map,newOrder, INFLOPlusNNdist, data_n_rows);
//            stratifier(INFLOVector,NNdistVector,INFLOPlusNNdist,map,newOrder,layerVector,data_n_rows);
//
//            Rcpp::IntegerVector clusteringRes = calc_cluster(calIsk,data_n_rows,k,nOfIsk,border,map);
//            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes, Rcpp::Named("layer") = layerVector, Rcpp::Named("border") = border);
//        }
//    }
}