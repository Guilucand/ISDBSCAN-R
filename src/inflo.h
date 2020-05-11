//
// Created by andrea on 14/04/20.
//

#ifndef ISDBSCAN_R_MASTER_INFLO_H
#define ISDBSCAN_R_MASTER_INFLO_H

#include "structs.h"

std::vector<double> calculate_inflo(std::vector<knn_dists<double>> const& knn_dists, std::vector<knn_influence_space> const& isk, double &max_inflo);

void calculate_ainflo(std::vector<int> const& nof_rnn, std::vector<knn_influence_space> const& isk, std::vector<double> &inflo);

void normalize_ainflo(std::vector<double> &inflo, double max_ainflo, double max_nndist);

#endif //ISDBSCAN_R_MASTER_INFLO_H
