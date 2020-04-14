//
// Created by andrea on 29/03/20.
//

#include "influence_space.h"
#include <iostream>
#include <algorithm>

// Calcolo i nodi influence_space (old ISK)
std::vector<knn_influence_space>
find_influence_space(std::vector<knn_dists<double>> const &dists, int k) {

    std::cout << "Finding influence_space components... " << std::endl;
    std::vector<knn_influence_space> influence_space(dists.size());

    for (int i = 0; i < dists.size(); i++) {
        knn_influence_space bic;
        for (int j = 0; j < k; j++) {
            auto el_index = dists[i][j].index;
            auto match = std::find_if(dists[el_index].begin(), dists[el_index].end(),
                                      [i](auto const &el) { return el.index == i; });

            if (match != dists[el_index].end()) {
                bic.push_back(el_index);
            }
        }
        influence_space[i] = bic;
    }
    std::cout << "Done." << std::endl;
    return influence_space;
}
