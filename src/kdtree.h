#ifndef ISDBSCAN_R_MASTER_KDTREE_H
#define ISDBSCAN_R_MASTER_KDTREE_H

#include <cmath>
#include "nearest_queue.h"

#define BRUTEFORCE_SIZE 64

template<typename T>
class kd_box {

    kd_point<T> mins;
    kd_point<T> maxs;
    int dimension;

    std::vector<kd_box<T>> childs;
    std::vector<kd_point<T>> points;

    template <typename It>
    double _compute_stdev(It begin, It end, int col) {
        It it = begin;
        double sum = 0.0;
        double sumsq = 0.0;
        double elcount = 0.0;
        while (it != end) {
            T val = it->features[col];
            sum += val;
            sumsq += val * val;
            ++elcount;
            ++it;
        }
        return sqrt((sumsq - sum) / (elcount - 1.0));
    }

public:
    template <typename It>
    kd_box(It begin, It end) {

        this->dimension = -1;

        size_t len = end - begin;
        if (len == 0) return;

        size_t clen = (*begin).features.size();
        this->mins.features.resize(clen, std::numeric_limits<T>::max());
        this->maxs.features.resize(clen, std::numeric_limits<T>::min());

        if (len < BRUTEFORCE_SIZE) {
            this->points.insert(this->points.begin(), begin, end);

            It it = begin;
            while (it != end) {
                for (int i = 0; i < clen; i++) {
                    this->mins.features[i] = std::min(this->mins.features[i], it->features[i]);
                    this->maxs.features[i] = std::max(this->maxs.features[i], it->features[i]);
                }
                ++it;
            }
            return;
        }

        int dimcnt = (*begin).features.size();

        std::pair<double, int> stdev = {0, -1};

        for (int i = 0; i < dimcnt; i++) {
            stdev = max(stdev, {_compute_stdev(begin, end, i), i});
        }


        int cmax = stdev.second;
        this->dimension = cmax;
        sort(begin, end, [cmax](auto &a, auto &b) {
            return a.features[cmax] < b.features[cmax];
        });

        int boxes_count = std::max(2, (int)ceil(log2(len) / log2(8)));

        for (int i = 0; i < boxes_count; i++) {
            It pstart = begin + len * i / boxes_count;
            It pend = begin + len * (i + 1) / boxes_count;
            this->childs.emplace_back(pstart, pend);
        }
        for (auto &child : childs) {
            for (int i = 0; i < clen; i++) {
                this->mins.features[i] = std::min(this->mins.features[i], child.mins.features[i]);
                this->maxs.features[i] = std::max(this->maxs.features[i], child.maxs.features[i]);
            }
        }
    }

    double distance_fast(kd_point<T> const& pt) const {
        if (this->dimension < 0) return _distance_slow(pt);
        T val = 0;
        val += std::max((T)0, mins.features[this->dimension] - pt.features[this->dimension]);
        val += std::max((T)0, pt.features[this->dimension] - maxs.features[this->dimension]);
        return val * val;
    }

    double _distance_slow(kd_point<T> const& pt) const {
        T dist = 0;
        for (int i = 0; i < pt.features.size(); i++) {
            T val = 0;
            val += std::max((T)0, mins.features[i] - pt.features[i]);
            val += std::max((T)0, pt.features[i] - maxs.features[i]);
            dist += val * val;
        }
        return dist;
    }

    bool distance_less_than(kd_point<T> const& pt, double threshold) const {

        if (distance_fast(pt) > threshold) return false;

        T dist = 0;
        for (int i = 0; i < pt.features.size(); i++) {
            T val = 0;
            val += std::max((T)0, mins.features[i] - pt.features[i]);
            val += std::max((T)0, pt.features[i] - maxs.features[i]);
            dist += val * val;
            if (dist > threshold) return false;
        }
        return true;
    }

    void search(kd_nearest<T>& dists) const {

        static int cnt = 0;
        cnt++;

        for (auto &pt : points) {
            dists.add_element(pt);
        }

        T _max = dists.get_maximum();

        std::vector<const kd_box<T>*> boxes; boxes.resize(childs.size());
        for (int i = 0; i < boxes.size(); i++) {
            boxes[i] = &childs[i];
        }

        sort(boxes.begin(), boxes.end(), [&](auto const& a, auto const& b) {
            return a->distance_fast(dists.x) < b->distance_fast(dists.x);
        });

        for (auto const &box : boxes) {
            if (box->distance_less_than(dists.x, _max)) {
                box->search(dists);
                _max = dists.get_maximum();
            }
        }
    }
};

#endif //ISDBSCAN_R_MASTER_KDTREE_H
