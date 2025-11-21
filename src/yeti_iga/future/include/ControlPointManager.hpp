#pragma once
#include <vector>
#include <cstddef>
#include <mutex>

struct ControlPointManager {
    int dim_phys = 3;

    // Contiguous storage : [x0, y0, z0, x1, y1, z1, ...]
    std::vector<double> coords;


    std::vector<double> field; // a field defined at control points

    // minimal thread-safety
    std::mutex mtx;

    size_t n_points() const {return coords.size() / dim_phys; }

    // reserve / append / set
    size_t add_point(const std::vector<double>& p) {
        std::lock_guard<std::mutex> lk(mtx);
        size_t id = n_points();
        coords.insert(coords.end(), p.begin(), p.end());
        return id;
    }

    // fill contiguous block from python pointer (fast import)
    void reserve_points(size_t n, int dim) {
        std::lock_guard<std::mutex> lk(mtx);
        dim_phys = dim;
        coords.reserve(n * dim);
    }
};