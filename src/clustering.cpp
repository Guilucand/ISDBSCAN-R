//
// Created by andrea on 29/03/20.
//

#include "clustering.h"
#include <iostream>

struct expand_cluster_stack {
    int pos_init;
    knn_influence_space::const_iterator it;
};

// Expand cluster iterative (to avoid stack overflows)
static void expand_cluster_iter(std::vector<knn_influence_space> const& kneighbours, std::vector<int> &all_clusters,
                       std::vector<int> &current_cluster, int pos_init, int cluster, int k,
                       int &assigned, std::vector<bool> &border) {

    std::vector<expand_cluster_stack> stack;
    stack.push_back({pos_init, kneighbours[pos_init].begin()});

    while (!stack.empty()) {

        bool do_recurse = false;

        expand_cluster_stack &frame = stack.back();

        if (kneighbours[frame.pos_init].size() > (k * 2.0 / 3.0)) {
            while (frame.it != kneighbours[frame.pos_init].end()) {
                auto kidx = *frame.it;
                if (all_clusters[kidx] == -1) {
                    current_cluster.push_back(kidx);
                    all_clusters[kidx] = cluster;
                    assigned++;

                    ++frame.it;
                    // Do recursion call
                    stack.push_back({kidx, kneighbours[kidx].begin()});
                    do_recurse = true;
                    break;
                }
                else {
                    ++frame.it;
                }
            }
        } else {
            border[frame.pos_init] = true;
        }
        if (!do_recurse) {
            stack.pop_back();
        }
    }
}

// Expand cluster
static void expand_cluster_rec(std::vector<knn_influence_space> const& kneighbours, std::vector<int> &all_clusters,
                               std::vector<int> &current_cluster, int pos_init, int cluster, int k,
                               int &assigned, std::vector<bool> &border) {
    if (kneighbours[pos_init].size() > (k * 2.0 / 3.0)) {
        for (auto &iskIdx : kneighbours[pos_init]) {
            if (all_clusters[iskIdx] == -1) {
                current_cluster.push_back(iskIdx);
                all_clusters[iskIdx] = cluster;
                assigned--;
                expand_cluster_rec(kneighbours, all_clusters, current_cluster, iskIdx, cluster, k,
                                   assigned, border);
            }
        }
    } else {
        border[pos_init] = true;
    }
}

std::vector<int>
calc_cluster(std::vector<knn_influence_space> const& kneighbours, int k, std::vector<bool> &border, std::vector<int> const& ranks_map, std::vector<bool> noise) {

    std::cout << "Computing clusters..." << std::endl;
    int nsamples = kneighbours.size();

    std::vector<int> all_clusters(nsamples, -1);

    // Set as noise
    for (int i = 0; i < noise.size(); i++) {
        if (noise[i]) all_clusters[i] = -2;
    }

    int cluster = 0;

    int assigned = 0;

    for (int idx_pos_init = 0; idx_pos_init < nsamples; idx_pos_init++) {

        int mapped_pos = ranks_map[idx_pos_init];

        if (all_clusters[mapped_pos] != -1) continue;
        all_clusters[mapped_pos] = cluster;
        assigned++;

        int just_assigned = 0;
        std::vector<int> current_cluster;
        current_cluster.push_back(mapped_pos);

//        expand_cluster_rec(kneighbours, all_clusters, current_cluster, pos_init, cluster, k, just_assigned, border);
        expand_cluster_iter(kneighbours, all_clusters, current_cluster, mapped_pos, cluster, k, just_assigned, border);

        if (just_assigned >= k) {
            std::cout << "Assigned cluster " << cluster << " with " << just_assigned << " points" << std::endl;
            assigned += just_assigned;
            cluster++;
        }
        else {
            for (auto &idxPoint : current_cluster) {
                if (idxPoint != mapped_pos) {
                    if (kneighbours[idxPoint].empty()) {
                        std::cout << "Error!" << std::endl;
                        all_clusters[idxPoint] = -1;
                        assigned--;
                    }
                }
                // TODO: Check if it's compliant with the algorithm
                all_clusters[idxPoint] = -1;
            }
            if (just_assigned > 0) {
                std::cout << "Not assigned cluster " << cluster << " with "
                          << just_assigned << " points" << std::endl;
            }
            all_clusters[mapped_pos] = -2;
        }
    }

    std::cout << "Done." << std::endl;
    return all_clusters;
}
