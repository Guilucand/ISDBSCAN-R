//
// Created by andrea on 29/03/20.
//

#ifndef ISDBSCAN_R_MASTER_CLUSTERING_H
#define ISDBSCAN_R_MASTER_CLUSTERING_H

#include "structs.h"

std::vector<int>
    calc_cluster(std::vector<knn_influence_space> const& kneighbours, int k, std::vector<bool> &border);

#endif //ISDBSCAN_R_MASTER_CLUSTERING_H
