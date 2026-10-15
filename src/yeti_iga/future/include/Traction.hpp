#pragma once
#include <Eigen/Dense>

// Value of a distributed boundary load (force per unit length along the
// edge), evaluated at a physical point. ConstantTraction is the only
// concrete kernel for now -- evaluate() already takes the physical point
// (not just internal parameters) so that a future Python-callback-based
// subclass (mirroring LocalOperator's pybind11 trampoline) can be added
// without changing PatchIntegrator's boundary-load integration loop.
class Traction {
public:
    virtual ~Traction() = default;

    virtual Eigen::Vector2d evaluate(const Eigen::Vector2d& physical_point) const = 0;
};

class ConstantTraction : public Traction {
public:
    explicit ConstantTraction(const Eigen::Vector2d& value) : value_(value) {}

    Eigen::Vector2d evaluate(const Eigen::Vector2d& /*physical_point*/) const override {
        return value_;
    }

private:
    Eigen::Vector2d value_;
};

// 3D equivalent of Traction/ConstantTraction, for a boundary load on a
// 3D solid patch's FACE (not an edge -- see
// PatchIntegrator::integrateBoundaryLoad3D()). A genuinely separate type
// from Traction (not a template/variant): this codebase keeps the 2D and 3D
// solid code paths textually separate throughout (see PatchIntegrator.hpp).
class Traction3D {
public:
    virtual ~Traction3D() = default;

    virtual Eigen::Vector3d evaluate(const Eigen::Vector3d& physical_point) const = 0;
};

class ConstantTraction3D : public Traction3D {
public:
    explicit ConstantTraction3D(const Eigen::Vector3d& value) : value_(value) {}

    Eigen::Vector3d evaluate(const Eigen::Vector3d& /*physical_point*/) const override {
        return value_;
    }

private:
    Eigen::Vector3d value_;
};
