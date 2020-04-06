//
// Created by andrea on 29/03/20.
//

#ifndef ISDBSCAN_R_MASTER_KNN_H
#define ISDBSCAN_R_MASTER_KNN_H

#include "structs.h"

std::vector<knn_dists<double>>
calc_knn_dists_approx(std::vector<kd_point<double>> const& data, int k, int trees, int checks);

std::vector<knn_dists<double>>
calc_knn_dists_exact_lowdim(std::vector<kd_point<double>> const& data, int k);

std::vector<knn_dists<double>>
calc_knn_dists_exact_highdim(std::vector<kd_point<double>> const& data, int k);

inline std::vector<knn_dists<double>>
calc_knn_dists(std::vector<kd_point<double>> const& data, int k, bool approx) {
    if (approx) {
        return calc_knn_dists_approx(data, k, 4, 128);
    } else {
        if (k < 50)
            return calc_knn_dists_exact_lowdim(data, k);
        else
            return calc_knn_dists_exact_highdim(data, k);
    }
}

#endif //ISDBSCAN_R_MASTER_KNN_H
