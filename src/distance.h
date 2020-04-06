#ifndef ISDBSCAN_R_MASTER_DISTANCE_H
#define ISDBSCAN_R_MASTER_DISTANCE_H

#include "structs.h"
#include <cassert>

template <typename T>
T get_distance(kd_point<T> const& a, kd_point<T> const& b) {
    assert(a.features.size() == b.features.size());
    int size = a.features.size();
    T dist = 0;
    for (int i = 0; i < size; i++) {
        T idist = a.features[i] - b.features[i];
        dist += idist * idist;
    }
    return dist;
}

#endif //ISDBSCAN_R_MASTER_DISTANCE_H
