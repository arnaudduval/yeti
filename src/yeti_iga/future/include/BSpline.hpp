#pragma once
#include <vector>
#include <pybind11/numpy.h>


namespace py = pybind11;

class BSpline {
public:
    BSpline(int degree_, py::array_t<double>kvector_);
    py::array_t<double> kvView() const;
    int FindSpan(double u) const;
    py::array_t<double> BasisFuns(int span, double u) const;
    double OneBasisFun(double u, int i);
    int getDegree() const {return degree; }


private:
    int degree;
    std::vector<double> kvector;
};
