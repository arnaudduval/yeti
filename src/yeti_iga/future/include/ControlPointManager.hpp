#pragma once
#include <algorithm>
#include <vector>
#include <cstddef>
#include <mutex>
#include <stdexcept>

struct ControlPointManager {
    int dim_phys = 3;

    // Contiguous storage : [x0, y0, z0, x1, y1, z1, ...]
    std::vector<double> coords;

    std::vector<double> field; // a field defined at control points

    // NURBS rational weights, indexed by global control point id.
    // Empty when all weights are 1.0 (pure B-spline, zero overhead).
    // Populated as soon as any weight != 1.0 is set; size always equals
    // n_points() when non-empty.
    std::vector<double> weights;

    // minimal thread-safety
    std::mutex mtx;

    size_t n_points() const {return coords.size() / dim_phys; }

    bool is_rational() const { return !weights.empty(); }

    // reserve / append / set
    // w: NURBS weight (default 1.0 = pure B-spline, no cost).
    // When w != 1.0, or when is_rational() is already true, the weights
    // vector is activated/extended so that every point has an explicit
    // weight (existing points without a set weight default to 1.0).
    size_t add_point(const std::vector<double>& p, double w = 1.0) {
        std::lock_guard<std::mutex> lk(mtx);
        size_t id = n_points();
        coords.insert(coords.end(), p.begin(), p.end());
        if (w != 1.0 || !weights.empty()) {
            // activate rational mode: fill any missing previous weights with 1.0
            if (weights.size() < id) weights.resize(id, 1.0);
            weights.push_back(w);
        }
        return id;
    }

    // Set the weight of an existing control point.
    // Activates rational mode (initialises every previous weight to 1.0 if
    // needed). Only useful post-construction; prefer add_point(p, w) when
    // building new geometries.
    void set_weight(size_t id, double w) {
        std::lock_guard<std::mutex> lk(mtx);
        size_t n = n_points();
        if (id >= n)
            throw std::out_of_range("ControlPointManager::set_weight: id out of range.");
        if (weights.empty()) weights.assign(n, 1.0);
        weights[id] = w;
        if (w == 1.0) {
            // check if all weights are now 1.0 again -> revert to B-spline fast path
            bool all_one = std::all_of(weights.begin(), weights.end(),
                                       [](double wi) { return wi == 1.0; });
            if (all_one) weights.clear();
        }
    }

    // Return the weight of a control point (1.0 if not rational, no lookup).
    double get_weight(size_t id) const {
        if (weights.empty()) return 1.0;
        return weights[id];
    }

    // fill contiguous block from python pointer (fast import)
    void reserve_points(size_t n, int dim) {
        std::lock_guard<std::mutex> lk(mtx);
        dim_phys = dim;
        coords.reserve(n * dim);
    }
};
