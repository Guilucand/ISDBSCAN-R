#include "knn.h"
#include "kdtree.h"
#include <iostream>

std::vector<knn_dists<double>>
        calc_knn_dists_exact_lowdim(std::vector<kd_point<double>> const& data, int k) {

    std::cout << "Computing knn distances lowdim..." << std::endl;

    std::vector<kd_point<double>> dataset_cpy = data;

    kd_box<double> box = kd_box<double>(dataset_cpy.begin(), dataset_cpy.end());

    size_t nsamples = data.size();
    std::vector<knn_dists<double>> result(nsamples);

    for (int i = 0; i < nsamples; i++) {
        kd_nearest<double> dists = kd_nearest<double>(data[i], k);
        box.search(dists);

        result[i] = dists.get_distances();

        if (i % 1000 == 0) {
            std::cout << "Knn finder... " << i << "/" << nsamples << std::endl;
        }
    }
    std::cout << "Knn finder... " << nsamples << "/" << nsamples << std::endl;
    std::cout << "Done." << std::endl;
    return result;
}

std::vector<knn_dists<double>>
        calc_knn_dists_exact_highdim(std::vector<kd_point<double>> const& data, int k) {

    std::cout << "Computing knn distances highdim..." << std::endl;

    std::vector<kd_point<double>> dataset_cpy = data;

    size_t nsamples = data.size();
    std::vector<knn_dists<double>> result(nsamples);

    for (int i = 0; i < nsamples; i++) {
        kd_nearest<double> dists = kd_nearest<double>(data[i], k);

        for (int j = 0; j < nsamples; j++) {
            dists.add_element(data[j]);
        }

        result[i] = dists.get_distances();

//        for (auto x : result[i]) {
//            std::cout << "(" << x.distance << "->" << x.index << ") ";
//        }
//        std::cout << std::endl;

        if (i % 1000 == 0) {
            std::cout << "Knn finder... " << i << "/" << nsamples << std::endl;
        }
    }
    std::cout << "Knn finder... " << nsamples << "/" << nsamples << std::endl;
    std::cout << "Done." << std::endl;
    return result;
}
