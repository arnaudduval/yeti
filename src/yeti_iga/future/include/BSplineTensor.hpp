#pragma once
#include "BSpline.hpp"
#include <pybind11/stl.h>

class BSplineTensor {
public:
    std::vector<BSpline> components;

    BSplineTensor(const std::vector<BSpline>& comps) : components(comps) {}
    py::array_t<int> FindSpanND(const py::array_t<double>& u) const;
    py::array_t<double> BasisFunsND(const py::array_t<int>& span, const py::array_t<double>& u) const;

    // Optimized versions
    void FindSpanND(const double* u, int* spans) const;
    // py::array_t<double> BasisFunsND(const int* spans, const double* u) const;
    

};

class BSplineSurface : public BSplineTensor {
public:
    BSplineSurface(const BSpline& su,
                   const BSpline& sv) : BSplineTensor({su, sv}) {}
};

class BSplineVolume : public BSplineTensor {
public:
    BSplineVolume(const BSpline& su,
                  const BSpline& sv,
                  const BSpline& sw) : BSplineTensor({su, sv, sw}) {}
};
