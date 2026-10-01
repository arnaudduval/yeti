#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BSpline.hpp"

// Weighted-quadrature (WQ) data for one parametric direction: a reduced set
// of quadrature points (fewer than standard Gauss for the same target
// accuracy), each basis function carrying its OWN weight at each point
// (unlike Gauss, where every basis function shares the same per-point
// weight) -- precomputed once by solving a small constrained least-squares
// problem per basis function row.
//
// Direct C++ port of pymfiga's
// src/yeti_iga/pymfiga/common/numerics/quadrature_rules/weighted_quadrature.py,
// both variants ("method 1" and "method 2"): they share the same point-
// placement scheme (midpoint rule, method-dependent point counts) and the
// same per-row least-squares weight computation, differing only in the
// target space used as the least-squares constraint and in how many
// distinct weight matrices come out of it:
//   - method "2" (the only variant exercised by pymfiga's own validated
//     regression benchmarks, benchs/pymfiga/*/test_*.py, QUADCLASS="wq",
//     QUADTYPE="2"): target space = SAME degree as the original spline, with
//     interior knot multiplicity increased by 1
//     (Operations.increase_multiplicity_to_knotvector in pymfiga). Only 2
//     distinct least-squares solves per basis function; W00==W01, W10==W11.
//   - method "1": target space = degree-1, knot vector = the original knot
//     vector stripped of its first and last entry. 4 distinct least-squares
//     solves per basis function (W00, W01, W10, W11 all different): W00/W10
//     constrain against the ORIGINAL (test) space basis at the WQ points,
//     W01/W11 against the reduced target space.
//
// Reference: Calabro, Sangalli, Tani (2017), "Fast formation of isogeometric
// Galerkin matrices by weighted quadrature".
//
// Scope: non-periodic B-splines only (matching every existing `future`
// caller); the knot vector is assumed to span exactly [0, 1] (matching
// pymfiga's own assumption and every knot vector used elsewhere in `future`).
struct WeightedQuadrature1D {
    int degree = 0;
    int nbctrlpts = 0;
    std::string quadtype = "2";  // "1" or "2", matching pymfiga's own naming

    std::vector<double> quadpts;  // WQ point positions, parametric space

    // Global sparse basis matrices at the WQ points: B0(q,i) = N_i(u_q),
    // B1(q,i) = dN_i/du(u_q). Shape (nbquadpts x nbctrlpts).
    Eigen::SparseMatrix<double> B0, B1;

    // Global sparse WQ weight matrices, shape (nbctrlpts x nbquadpts): W_ab
    // is the weight of basis function i at quadrature point q used when
    // integrating a product of a test function differentiated `a` times
    // against a trial function differentiated `b` times (a,b in {0,1}).
    // W00: value*value (mass-like). W11: derivative*derivative
    // (stiffness-like). W01/W10: the mixed terms. For method "2" specifically
    // W00==W01 and W10==W11 (see above); for method "1" all four differ.
    Eigen::SparseMatrix<double> W00, W01, W10, W11;

    // Build the WQ rule for a given BSpline's degree/knot vector, using
    // pymfiga's method "1" or "2" (default "2", matching pymfiga's own
    // validated benchmarks).
    static WeightedQuadrature1D build(const BSpline& bspline, const std::string& quadtype = "2");
};
