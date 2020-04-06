//
// Created by andrea on 29/03/20.
//

#ifndef ISDBSCAN_R_MASTER_NEAREST_QUEUE_H
#define ISDBSCAN_R_MASTER_NEAREST_QUEUE_H

#include "structs.h"
#include "distance.h"
#include <queue>
#include <limits>
#include <algorithm>

template<typename T>
class kd_nearest {

    inline void _add_element(kd_point<T> const& elem, T dist) {
        points.push({dist, elem});
        min_distance = std::min(min_distance, dist);
        max_distance = points.top().distance;
    }
public:

    kd_nearest(kd_point<T> const& x, int kcount) {
        this->kcount = kcount+1;
        this->x = x;
    }

    struct kd_point_cmp {
        T distance;
        kd_point<T> pt;
        bool operator <(kd_point_cmp const& oth) const { return this->distance < oth.distance; }
    };

    int kcount;
    kd_point<T> x;
    T min_distance = std::numeric_limits<T>::max();
    T max_distance = std::numeric_limits<T>::min();
    std::priority_queue<kd_point_cmp> points;

    inline T get_maximum() {
        if (points.size() < kcount) {
            return std::numeric_limits<T>::max();
        }
        return max_distance;
    }

    inline void add_element(kd_point<T> const& elem) {
        T dist = get_distance(x, elem);
        if (points.size() < kcount) {
            _add_element(elem, dist);
            return;
        }
        else if (get_maximum() > dist) {
            points.pop();
            _add_element(elem, dist);
        }
    }

    std::vector<kd_point_cmp> get_elem() {
        auto pq = points;


        int sz = pq.size() - 1;
        std::vector<kd_point_cmp> result(std::max(0, sz));
        for (auto &pt : result) {
            pt = pq.top(); pq.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    knn_dists<T> get_distances() {
        auto pq = points;


        int sz = pq.size() - 1;
        knn_dists<double> result(std::max(0, sz));
        for (auto &pt : result) {
            auto top = pq.top(); pq.pop();
            pt = { top.distance, top.pt.index };
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

};


#endif //ISDBSCAN_R_MASTER_NEAREST_QUEUE_H
