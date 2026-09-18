#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Material.hpp"

// Kirchhoff-Love (rotation-free) shell kinematics: covariant/contravariant
// surface metric, curvature, and the isotropic membrane/bending constitutive
// matrix. Direct C++ translation of the legacy Fortran shell element
// (src/yeti_iga/main/shell/curvilinearCoordinates.f, UELMAT_shell.f,
// USFMEM_shell.f, USFBND_shell.f) -- see those files for the reference
// formulas this code reproduces.
//
// A Kirchhoff-Love shell has 3 translational DOFs per control point (no
// rotation field): all membrane and bending strain measures below are
// derived purely from first/second derivatives of the mid-surface position.

// Covariant surface quantities at one Gauss point, plus the derived
// contravariant metric. a1/a2 are the covariant tangents (first derivatives
// of the physical position w.r.t. the two parametric directions), a3 the
// unit normal. a11/a12/a22 are the second-derivative vectors needed for
// curvature (a12 = da1/dv = da2/du, by symmetry of mixed partials -- stored
// once, matching curvilinearCoordinates.f's dAI1dxi(2,:) == dAI2dxi(1,:)).
struct ShellGeometry {
    Eigen::Vector3d a1, a2, a3;
    Eigen::Vector3d a11, a12, a22;
    Eigen::Matrix2d AAI;   // covariant metric:    AAI(i,j) = a_i . a_j
    Eigen::Matrix2d AAE;   // contravariant metric: AAE = AAI^-1
    double Area;            // |a1 x a2| -- surface metric determinant sqrt(det(AAI))
};

// Builds a1, a2, a3, a11, a12, a22, AAI, AAE, Area from the physical-space
// gradients/second-derivatives of the basis functions and the 3D control
// points of the current span (u-fastest order, matching
// patch.control_points_for_span()). pts[a] must point at 3 contiguous
// physical doubles (x,y,z) for control point a.
ShellGeometry computeShellGeometry(
    const std::vector<const double*>& pts,
    const std::vector<double>& dRdu, const std::vector<double>& dRdv,
    const std::vector<double>& d2Rdu2, const std::vector<double>& d2Rdv2,
    const std::vector<double>& d2Rdudv);

// Membrane strain-displacement operator for one control point (3x3: rows =
// membrane strain components eps11, eps22, 2*eps12; columns = physical x,y,z
// DOF directions). Matches USFMEM_shell.f's BoJ (per node, transposed to
// build BoI): K_ab = membraneB(a)^T * matH_membrane * membraneB(b).
Eigen::Matrix3d membraneB(const ShellGeometry& g, double dRdu_a, double dRdv_a);

// Bending (curvature) strain-displacement operator for one control point
// (3x3: rows = curvature components kappa11, kappa22, 2*kappa12; columns =
// physical x,y,z). Matches USFBND_shell.f's BoJ.
Eigen::Matrix3d bendingB(const ShellGeometry& g, double dRdu_a, double dRdv_a,
                         double d2Rdu2_a, double d2Rdv2_a, double d2Rdudv_a);

// Isotropic Kirchhoff-Love shell constitutive law: builds the 3x3 "matH"
// material matrix from the contravariant metric at a Gauss point (matH
// varies over a curved shell, unlike a flat-plate material matrix). The
// caller scales matH by thickness (membrane) or thickness^3/12 (bending) --
// same matH, no membrane-bending coupling term, matching the legacy
// single-homogeneous-isotropic-layer assumption exactly
// (UELMAT_shell.f:150-173).
class KirchhoffLoveShellLaw {   // NOT a ConstitutiveLaw subclass -- see class
                                 // docstring in ConstitutiveLaw.hpp for why.
public:
    explicit KirchhoffLoveShellLaw(const Material& m) : mat_(m) {}
    const Material& material() const { return mat_; }

    Eigen::Matrix3d matH(const Eigen::Matrix2d& AAE) const {
        const double nu = mat_.nu;
        const double coef = mat_.E / (1.0 - nu * nu);
        Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
        H(0,0) = AAE(0,0) * AAE(0,0);
        H(1,1) = AAE(1,1) * AAE(1,1);
        H(2,2) = 0.5 * ((1.0 - nu) * AAE(0,0) * AAE(1,1)
                       + (1.0 + nu) * AAE(0,1) * AAE(0,1));
        H(0,1) = H(1,0) = nu * AAE(0,0) * AAE(1,1) + (1.0 - nu) * AAE(0,1) * AAE(0,1);
        H(0,2) = H(2,0) = AAE(0,0) * AAE(0,1);
        H(1,2) = H(2,1) = AAE(1,1) * AAE(0,1);
        return coef * H;
    }

private:
    Material mat_;
};
