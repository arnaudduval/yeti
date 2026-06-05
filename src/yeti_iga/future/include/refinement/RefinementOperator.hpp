#pragma once
#include <memory>
#include <Eigen/Dense>
#include "Patch.hpp"

class RefinementOperator {
public:
    virtual ~RefinementOperator() = default;

    // Refine a patch in-place.
    // transition_matrix : (nb_new_cp x nb_old_cp)
    virtual void refine(
        Patch& patch,
        Eigen::MatrixXd& transition_matrix
    ) const = 0;

    virtual std::string getType() const = 0;
};
