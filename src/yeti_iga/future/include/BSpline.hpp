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
    void BasisFuns_raw(int span, double u, double* out) const;
    // TODO plutôt appeler cette fonction ..._raw ?
    void BasisFunsDerivatives(int span, double u, int n, double* ders_out) const;
    double OneBasisFun(double u, int i);
    int getDegree() const {return degree; }
    const std::vector<double>& getKnotVector() const { return kvector; }


private:
    int degree;
    std::vector<double> kvector;
};

// Interface for python binding
py::array_t<double> bspline_basis_funs_derivatives(const BSpline& self, int span, double u, int d);
