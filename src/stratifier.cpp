//
// Created by andrea on 14/04/20.
//

#include <numeric>
#include <algorithm>
#include "stratifier.h"

bool create_partition(
        std::vector<double> const& inflo,
        std::vector<std::pair<double, int>> const& inflo_plus_nn_sorted,
        std::vector<int> const& ranks_map,
        std::vector<int> &layer_vec,
        const double mean_ipnns,
        int &cut,
        const int layer,
        const double inflo_threshold){

    int new_cut = 0;

    for(int i = 0;  i < inflo.size(); i++ ){
        if(inflo_plus_nn_sorted[i].first < mean_ipnns){
            new_cut = i;
        }
    }

    double interval_mean_inflo = 0;

    for (int i = cut; i < new_cut + 1; i++ ){
        interval_mean_inflo += inflo[ranks_map[i]];
    }
    interval_mean_inflo = interval_mean_inflo / (new_cut - cut + 1);

    for (int i = cut; i < new_cut +1; i++){
        layer_vec[ranks_map[i]] = layer;
    }

    cut = new_cut + 1;
    return interval_mean_inflo >= inflo_threshold;
}


// Stratifier
std::vector<int> stratifier(std::vector<double> const& inflo,
        std::vector<std::pair<double, int>> const& inflo_plus_nn_sorted,
        std::vector<int> const& ranks_map) {

    std::vector<int> layer_vec(inflo.size());

    double sum_inflo = std::accumulate(inflo.begin(), inflo.end(), 0.0);
    double mean_inflo = sum_inflo / inflo.size();

    std::vector<double> diff(inflo.size());
    std::transform(inflo.begin(), inflo.end(), diff.begin(), [mean_inflo](double x) { return x - mean_inflo; });

    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double var_inflo = (sq_sum / inflo.size());

    int layer =  0;
    int cut = 0;

    std::vector<int> cuts;
    double inflo_treshold = mean_inflo + var_inflo;

    //begin partitioning
    while(cut < inflo.size()) {
        double mean_ipnns = 0;

        for (int i = cut; i < inflo.size(); i++) {
            mean_ipnns += inflo_plus_nn_sorted[i].first;
        }
        mean_ipnns /= (double)(inflo.size() - cut);

        bool stop = create_partition(inflo, inflo_plus_nn_sorted, ranks_map, layer_vec, mean_ipnns, cut, layer, inflo_treshold);

        layer++;
        cuts.push_back(cut);

        if(stop) {

            for(int i = cut; i < inflo.size(); i++) {
                layer_vec[ranks_map[i]] = layer;
            }

            cuts.push_back(cut);
            break;
        }
    }
    return layer_vec;
}
