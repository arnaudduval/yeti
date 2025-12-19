#pragma once
#include <vector>
#include <cstddef>
#include "IGABasis1D.hpp"
#include "Patch.hpp"
#include <Eigen/Dense>


struct MaterialProperties {
    double E;
    double nu;
    double thickness;   // Thickness (for plane problems)
};

struct ElementMatrix {
    // global indices for the local basis
    std::vector<size_t> global_indices;
    // dense local stiffness
    Eigen::MatrixXd K;
    size_t nb_loc = 0;

    // Constructor (initialize K with required shape)
    ElementMatrix(size_t size) : nb_loc(size), K(2*size, 2*size) {
        K.setZero();
    }
};

class IGAAssembler2D {
public:
    // Construct with the patch and pre-computed per-direction IGABasis1D
    IGAAssembler2D(const Patch& patch,
                   const IGABasis1D& basis_u,
                   const IGABasis1D& basis_v)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v)
    {
        // sanity checks
        assert(patch_.local_shape.size() == 2);
        p_u_ = patch_.tensor.components[0].getDegree();
        p_v_ = patch_.tensor.components[1].getDegree();
        nb_loc_ = (p_u_+1) * (p_v_+1);

        // Hard coding of material properties, TODO : should be set properly
        material_.E = 210000.;
        material_.nu = 0.3;
    }

    // Assemble per-element stiffness blocks and return them
    // (it DOES NOT assemble into a global sparse matrix)
    std::vector<ElementMatrix> assemble_stiffness() const;

    Eigen::Matrix3d computeConstitutiveMatrix() const;

private:
    const Patch& patch_;
    const IGABasis1D& basis_u_;
    const IGABasis1D& basis_v_;
    int p_u_, p_v_;
    size_t nb_loc_;
    MaterialProperties material_;       // Should be defined as a member of patch ?
};
