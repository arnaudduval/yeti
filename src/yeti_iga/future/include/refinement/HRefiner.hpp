#pragma once
#include "refinement/RefinementOperator.hpp"


class HRefiner : public RefinementOperator {
public:
    ~HRefiner() override = default;
    HRefiner(int direction, double knot) : direction_(direction), knot_(knot) {}

    // Refine patch in-place by inserting one knot (h-refinement / Boehm's algorithm).
    // Only the p blended CPs per line are added to cp_manager — unchanged CPs are reused.
    void refine(
        Patch& patch,
        Eigen::MatrixXd& transition_matrix
    ) const override;

    std::string getType() const override { return "h"; }

private:
    int direction_;
    double knot_;
};
