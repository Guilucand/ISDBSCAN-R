//
// Created by andrea on 14/04/20.
//

#ifndef ISDBSCAN_R_MASTER_STRATIFIER_H
#define ISDBSCAN_R_MASTER_STRATIFIER_H

#include "structs.h"

std::vector<int> stratifier(std::vector<double> const& inflo,
                            std::vector<std::pair<double, int>> const& inflo_plus_nn_sorted,
                            std::vector<int> const& ranks_map);

#endif //ISDBSCAN_R_MASTER_STRATIFIER_H
