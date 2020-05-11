#include "knn.h"
#include <iostream>
#include <flann/flann.h>

std::vector<knn_dists<double>>
        calc_knn_dists_approx(std::vector<kd_point<double>> const& data, int k, int trees, int checks) {

    std::cout << "Computing approximate knn distances..." << std::endl;

    std::vector<double> flat_dataset;
    flat_dataset.reserve(data.size() * data[0].features.size());
    for (auto &row : data) {
        flat_dataset.insert(flat_dataset.end(), row.features.begin(), row.features.end());
    }

    auto flann_dataset = flann::Matrix<double>(flat_dataset.data(), data.size(),
                                               data[0].features.size());

    flann::Index <flann::L2<double>> index(flann_dataset, flann::KDTreeIndexParams(trees));

    std::cout << "Building flann index..." << std::endl;
    index.buildIndex();


    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<double>> dists_flann;


    std::cout << "Computing distances..." << std::endl;
    index.knnSearch(flann_dataset, indices, dists_flann, k + 1, flann::SearchParams(checks));

    size_t nsamples = data.size();
//    size_t nfeatures = !data.empty() ? data.front().features.size() : 0;

    std::vector<knn_dists<double>> result(nsamples);

    for (size_t i = 0; i < nsamples; i++) {
        result[i].reserve(k);
        for (int j = 0; j < k; j++) {
            result[i].push_back(knn_dist<double>{
                    dists_flann[i][j + 1],
                    indices[i][j + 1],
            });
        }
    }

    std::cout << "Done." << std::endl;
    return result;
}