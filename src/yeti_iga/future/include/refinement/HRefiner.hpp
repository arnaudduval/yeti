#pragma once
#include "refinement/RefinementOperator.hpp"


class HRefiner : public RefinementOperator {
public:
    ~HRefiner() override = default;
    HRefiner(int direction, double knot) : direction_(direction), knot_(knot) {}

    // Refine a single patch by inserting one knot (h-refinement / Boehm's algorithm)
    std::shared_ptr<Patch> refine(
        const Patch& patch,
        Eigen::MatrixXd& transition_matrix
    ) const override;

    std::string getType() const override { return "h"; }

private:
    int direction_;
    double knot_;

    void refine1DLine(
        const Patch& patch,
        int k,
        size_t line_start_old,
        size_t line_stride_old,
        size_t line_start_new,
        size_t line_stride_new,
        std::vector<std::vector<double>>& new_control_points,
        Eigen::MatrixXd& transition_matrix
    ) const;
};
