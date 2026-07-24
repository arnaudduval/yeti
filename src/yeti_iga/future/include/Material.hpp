#pragma once

// Physical material constants and derived elastic moduli.
// Immutable once constructed; all derived quantities are computed on demand.
struct Material {
    double E;
    double nu;
    double rho       = 0.0;   // mass density
    double thickness = 1.0;   // out-of-plane thickness for 2D plane problems

    double mu()           const { return E / (2.0 * (1.0 + nu)); }
    double lambda_3d()    const { return nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu)); }
    double bulk_modulus() const { return E / (3.0 * (1.0 - 2.0 * nu)); }

    static Material steel()     { return {210000.0, 0.30,  7850.0}; }
    static Material aluminium() { return {70000.0,  0.33,  2700.0}; }
    static Material concrete()  { return {30000.0,  0.20,  2400.0}; }
};
