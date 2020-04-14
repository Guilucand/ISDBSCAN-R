//
// Created by andrea on 29/03/20.
//

#ifndef ISDBSCAN_R_MASTER_INFLUENCE_SPACE_H
#define ISDBSCAN_R_MASTER_INFLUENCE_SPACE_H

#include "structs.h"

// Calcolo i nodi influence_space (old ISK)
std::vector<knn_influence_space>
        find_influence_space(std::vector<knn_dists<double>> const &dists, int k);

#endif //ISDBSCAN_R_MASTER_INFLUENCE_SPACE_H
