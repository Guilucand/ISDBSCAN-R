//
// Created by andrea on 29/03/20.
//

#include "biconnected.h"
#include <iostream>
#include <algorithm>

// Calcolo i nodi biconnected (old ISK)
std::vector<knn_biconnected>
find_biconnected(std::vector<knn_dists<double>> const &dists, int k) {

    std::cout << "Finding biconnected components... " << std::endl;
    std::vector<knn_biconnected> biconnected(dists.size());

    for (int i = 0; i < dists.size(); i++) {
        knn_biconnected bic;
        for (int j = 0; j < k; j++) {
            auto el_index = dists[i][j].index;
            auto match = std::find_if(dists[el_index].begin(), dists[el_index].end(),
                                      [i](auto const &el) { return el.index == i; });

            if (match != dists[el_index].end()) {
                bic.push_back(el_index);
            }
        }
        biconnected[i] = bic;
    }
    std::cout << "Done." << std::endl;
    return biconnected;
}
