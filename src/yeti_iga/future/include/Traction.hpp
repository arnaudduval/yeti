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
