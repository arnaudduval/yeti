#include "ConstitutiveLaw.hpp"

PhysMatrix IsotropicElastic::stiffness_density(
    const PhysVector& v,
    const PhysVector& w,
    const PhysVector& /*x_phys*/) const
{
    int n = static_cast<int>(v.size());
    return lambda_ * (v * w.transpose())
         + mu_     * (w * v.transpose())
         + (mu_ * v.dot(w)) * PhysMatrix::Identity(n, n);
}

PlaneStress::PlaneStress(const Material& m) : IsotropicElastic(m) {
    mu_     = m.mu();
    lambda_ = m.nu * m.E / (1.0 - m.nu * m.nu);
}

PlaneStrain::PlaneStrain(const Material& m) : IsotropicElastic(m) {
    mu_     = m.mu();
    lambda_ = m.lambda_3d();
}
