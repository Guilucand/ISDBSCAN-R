//
// Created by andrea on 14/04/20.
//

#ifndef ISDBSCAN_R_MASTER_STRATIFIER_H
#define ISDBSCAN_R_MASTER_STRATIFIER_H

#include "structs.h"

std::vector<int> stratifier(std::vector<double> const& ainflo,
                            std::vector<std::pair<double, int>> const& dfk_sorted,
                            std::vector<int> const& ranks_map,
                            int &noise_layer);

#endif //ISDBSCAN_R_MASTER_STRATIFIER_H
