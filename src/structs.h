#ifndef ISDBSCAN_R_MASTER_STRUCTS_H
#define ISDBSCAN_R_MASTER_STRUCTS_H

#include <vector>

template <typename T>
struct kd_point {
    std::vector<T> features;
    size_t index = -1;
};

template <class T>
struct knn_dist {
    T distance;
    size_t index;
    bool operator <(knn_dist<T> const& oth) const {
        return distance < oth.distance;
    }
};

using knn_influence_space = std::vector<int>;

template <class T>
using knn_dists = std::vector<knn_dist<T>>;

struct IsdbscanResult {
    std::vector<int> clusters;
    std::vector<int> layer;
    std::vector<bool> border;
};

#endif //ISDBSCAN_R_MASTER_STRUCTS_H
