//
// Created by andrea on 14/04/20.
//

#include "inflo.h"
#include "structs.h"

std::vector<double> calculate_inflo(std::vector<knn_dists<double>> const& knn_dists, std::vector<knn_influence_space> const& isk, double &max_inflo) {

    max_inflo = 0;
    std::vector<int> inf_inflo;
    std::vector<double> inflo_result(knn_dists.size());

    for (int i = 0; i < knn_dists.size(); i++){

        double current_knn_max = knn_dists[i].back().distance;

        int isk_num = isk[i].size();

        if (isk_num > 0) {
            double sum = 0.0;
            for (int j = 0; j < isk_num; j++){
                double knn_max = knn_dists[isk[i][j]].back().distance;
                double density_j = 1.0 / knn_max;
                sum += density_j;
            }

            double density_i = 1.0 / current_knn_max;
            sum = sum / (isk_num * density_i);
            inflo_result[i] = sum;
            max_inflo = std::max(max_inflo, sum);
        }
        else{
            inf_inflo.push_back(i);
        }
    }

    for (auto &idxPoint : inf_inflo) {
        inflo_result[idxPoint] = max_inflo;
    }

    return inflo_result;
}

void calculate_ainflo(std::vector<int> const& nof_rnn, std::vector<knn_influence_space> const& isk, std::vector<double> &inflo) {
    for (int i = 0; i < isk.size(); i++) {
        if (!isk[i].empty()) {
            inflo[i] /= nof_rnn[i];
        }
    }
}

// Calculate normalized AINFLO
void normalize_ainflo(std::vector<double> &inflo, double max_ainflo, double max_nndist) {
    for (int i = 0; i < inflo.size(); i++) {
        inflo[i] = (inflo[i] / max_ainflo) * max_nndist;
    }
}