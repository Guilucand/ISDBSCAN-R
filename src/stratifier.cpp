//
// Created by andrea on 14/04/20.
//

#include <numeric>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "stratifier.h"

int find_strata(
        std::vector<double> const& ainflo,
        std::vector<std::pair<double, int>> const& inflo_plus_nn_sorted,
        std::vector<int> const& ranks_map,
        std::vector<int> &layer_vec,
        const int cut,
        const int layer,
        const double inflo_threshold) {

    int new_cut = cut;

    double mean_df = 0;
    for (int i = cut; i < ainflo.size(); i++) {
        mean_df += inflo_plus_nn_sorted[i].first;
    }
    mean_df /= (double)(ainflo.size() - cut);


    while (inflo_plus_nn_sorted[new_cut].first < mean_df) {
        new_cut++;
    }

    double interval_mean_ainflo = 0;
    for (int i = cut; i < new_cut; i++ ){
        interval_mean_ainflo += ainflo[ranks_map[i]];
        layer_vec[ranks_map[i]] = layer;
    }
    interval_mean_ainflo = interval_mean_ainflo / (new_cut - cut);


    if (interval_mean_ainflo < inflo_threshold)
        return new_cut;

    std::cout << "Last layer!" << std::endl;
    for(int i = new_cut; i < ainflo.size(); i++) {
        layer_vec[ranks_map[i]] = layer;
    }
    return ainflo.size();
}


// Stratifier
std::vector<int> stratifier(std::vector<double> const& ainflo,
        std::vector<std::pair<double, int>> const& dfk_sorted,
        std::vector<int> const& ranks_map,
        int &noise_layer) {

    std::vector<int> layer_vec(ainflo.size());

    int layer =  0;
    int cut = 0;

    double sum_inflo = std::accumulate(ainflo.begin(), ainflo.end(), 0.0);
    double mean_inflo = sum_inflo / ainflo.size();

    std::vector<double> diff(ainflo.size());
    std::transform(ainflo.begin(), ainflo.end(), diff.begin(), [mean_inflo](double x) { return x - mean_inflo; });

    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double var_inflo = (sq_sum / ainflo.size());



    std::vector<int> cuts;
    double inflo_treshold = mean_inflo + var_inflo;

    //begin partitioning
    while(cut < ainflo.size() - 1) {

        cut = find_strata(ainflo, dfk_sorted, ranks_map, layer_vec, cut, layer,
                                inflo_treshold);
        layer++;
        cuts.push_back(cut);
    }
    noise_layer = layer - 1;

    std::vector<int> freqs(layer, 0);
    for (auto lay : layer_vec) {
        freqs[lay]++;
    }

    std::cout << "Layers distribution" << std::endl;
    for (int i = 0; i < freqs.size(); i++) {
        std::cout << "Number: " << i << "\t" << freqs[i] << std::endl;
    }

    return layer_vec;
}
