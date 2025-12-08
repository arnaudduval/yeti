#pragma once
#include <vector>
#include <cstddef>
#include "IGABasis1D.hpp"
#include "Patch.hpp"

struct ElementMatrix {
    // global indices for the local basis
    std::vector<size_t> global_indices;
    // dense local stiffness (nb_loc x nb_loc), row major
    std::vector<double> K;
    size_t nb_loc = 0;
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
    }

    // Assemble per-element stiffness blocks and return them
    // (it DOES NOT assemble into a global sparse matrix)
    std::vector<ElementMatrix> assemble_stiffness() const;

private:
    const Patch& patch_;
    const IGABasis1D& basis_u_;
    const IGABasis1D& basis_v_;
    int p_u_, p_v_;
    size_t nb_loc_;
};
