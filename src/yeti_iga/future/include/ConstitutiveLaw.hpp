#pragma once
#include <Eigen/Dense>
#include "Material.hpp"

// Abstract B-free constitutive law. Subclasses implement stiffness_density()
// using the identity C{v,w} = lambda*(v x w^T) + mu*(w x v^T) + mu*(v.w)*I
// (Planas, Romero & Sancho 2012, CMAME 217-220, 226-235).
// This avoids Voigt encoding, is dimension-agnostic, and is ~10x more efficient
// than B^T D B for isotropic materials (3 rank-1 updates vs one 3x3 product).
// Python subclassing is supported via the PyConstitutiveLaw trampoline in bindings.cpp.
class ConstitutiveLaw {
public:
    // Number of displacement DOFs per control point (2 for plane problems, 3 for 3D solid).
    virtual int n_dofs_per_cp() const = 0;

    // Elementary stiffness block K^{ab} = integral of stiffness_density(grad_a, grad_b, x_phys).
    // grad_a / grad_b: physical-space gradients of basis functions a and b (size = physical dim).
    // x_phys: physical coordinates of the Gauss point (needed for axisymmetric laws).
    // Returns a (n x n) matrix where n = n_dofs_per_cp().
    virtual Eigen::MatrixXd stiffness_density(
        const Eigen::VectorXd& grad_a,
        const Eigen::VectorXd& grad_b,
        const Eigen::VectorXd& x_phys) const = 0;

    const Material& material() const { return mat_; }
    virtual ~ConstitutiveLaw() = default;

protected:
    explicit ConstitutiveLaw(const Material& m) : mat_(m) {}
    Material mat_;
};

// Isotropic linear elastic laws: C{v,w} = lambda*(v x w^T) + mu*(w x v^T) + mu*(v.w)*I.
// PlaneStress and PlaneStrain differ only in lambda_ (mu_ is the same).
class IsotropicElastic : public ConstitutiveLaw {
protected:
    double lambda_, mu_;
    explicit IsotropicElastic(const Material& m) : ConstitutiveLaw(m) {}

public:
    Eigen::MatrixXd stiffness_density(
        const Eigen::VectorXd& v,
        const Eigen::VectorXd& w,
        const Eigen::VectorXd& /*x_phys*/) const override;
};

// 2D plane-stress: lambda_eff = nu*E / (1 - nu^2), n_dofs = 2.
class PlaneStress : public IsotropicElastic {
public:
    explicit PlaneStress(const Material& m);
    int n_dofs_per_cp() const override { return 2; }
};

// 2D plane-strain: lambda_eff = lambda_3D = nu*E / ((1+nu)(1-2nu)), n_dofs = 2.
class PlaneStrain : public IsotropicElastic {
public:
    explicit PlaneStrain(const Material& m);
    int n_dofs_per_cp() const override { return 2; }
};
