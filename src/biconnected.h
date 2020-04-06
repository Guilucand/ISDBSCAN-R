//
// Created by andrea on 29/03/20.
//

#ifndef ISDBSCAN_R_MASTER_BICONNECTED_H
#define ISDBSCAN_R_MASTER_BICONNECTED_H

#include "structs.h"

// Calcolo i nodi biconnected (old ISK)
std::vector<knn_biconnected>
        find_biconnected(std::vector<knn_dists<double>> const &dists, int k);

#endif //ISDBSCAN_R_MASTER_BICONNECTED_H
