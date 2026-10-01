#include "WQMatrixFreeStiffness.hpp"
#include "MatrixFree2D.hpp"
#include <stdexcept>

namespace {

// Extract the isotropic Lame parameters (lambda, mu) from any
// ConstitutiveLaw purely by numerically probing stiffness_density() with
// unit vectors -- no need for a new accessor on ConstitutiveLaw, and it
// works for PlaneStress/PlaneStrain (or any future IsotropicElastic
// subclass) alike. Derivation: for future's B-free identity
// K[I,J] = lambda*v[I]*w[J] + mu*w[I]*v[J] + mu*(v.w)*delta_IJ,
// evaluating at v=w=e0=(1,0) gives K[0,0]=lambda+2*mu, K[1,1]=mu.
std::pair<double, double> extract_lame(const ConstitutiveLaw& law) {
    PhysVector e0(2);
    e0 << 1.0, 0.0;
    PhysVector x(2);
    x << 0.0, 0.0;
    PhysMatrix K = law.stiffness_density(e0, e0, x);
    double mu = K(1, 1);
    double lambda = K(0, 0) - 2.0 * mu;
    return {lambda, mu};
}

} // namespace

WQMatrixFreeStiffness::WQMatrixFreeStiffness(
    const Patch& patch, const ConstitutiveLaw& law, const std::string& quadtype)
    : patch_(patch)
{
    if (law.n_dofs_per_cp() != 2)
        throw std::invalid_argument(
            "WQMatrixFreeStiffness: only 2D laws (n_dofs_per_cp() == 2, i.e. "
            "PlaneStress/PlaneStrain) are supported.");

    const BSpline& bsp_u = patch.tensor.components[0];
    const BSpline& bsp_v = patch.tensor.components[1];
    wq_u_ = WeightedQuadrature1D::build(bsp_u, quadtype);
    wq_v_ = WeightedQuadrature1D::build(bsp_v, quadtype);

    nbctrlpts_u_ = wq_u_.nbctrlpts;
    nbctrlpts_v_ = wq_v_.nbctrlpts;
    nq_u_ = static_cast<int>(wq_u_.quadpts.size());
    nq_v_ = static_cast<int>(wq_v_.quadpts.size());

    if (static_cast<int>(patch.local_shape[0]) != nbctrlpts_u_ ||
        static_cast<int>(patch.local_shape[1]) != nbctrlpts_v_)
        throw std::invalid_argument(
            "WQMatrixFreeStiffness: patch control-point grid size does not "
            "match the WeightedQuadrature1D rules built from its own "
            "BSplineTensor -- unexpected mismatch between patch.tensor and "
            "patch.local_shape.");
    if (!patch.dof_manager)
        throw std::invalid_argument("WQMatrixFreeStiffness: patch has no dof_manager.");

    const int nctrl = nbctrlpts_u_ * nbctrlpts_v_;
    const bool rational = patch.cp_manager->is_rational();

    Eigen::VectorXd x_ctrl(nctrl), y_ctrl(nctrl), w_ctrl;
    if (rational) w_ctrl.resize(nctrl);

    for (int jv = 0; jv < nbctrlpts_v_; ++jv) {
        for (int iu = 0; iu < nbctrlpts_u_; ++iu) {
            int idx = jv * nbctrlpts_u_ + iu;  // u-fastest, matches Patch::local_cp_ptr
            const double* P = patch.local_cp_ptr(static_cast<size_t>(idx));
            x_ctrl(idx) = P[0];
            y_ctrl(idx) = P[1];
            if (rational)
                w_ctrl(idx) = patch.cp_manager->get_weight(patch.global_indices[idx]);
        }
    }

    // Physical Jacobian components at every WQ point, gathered matrix-free
    // (Phase B) from control-point coordinates -- generalizes
    // PatchIntegrator::evaluateGaussPointGeometry(NURBS)'s per-span
    // quotient rule to the whole patch at once. Convention matches
    // PatchIntegrator exactly: J11=dx/du, J21=dy/du, J12=dx/dv, J22=dy/dv.
    Eigen::VectorXd J11, J12, J21, J22;
    if (!rational) {
        J11 = matrix_free_apply_2d(wq_u_.B1, wq_v_.B0, x_ctrl, false);
        J12 = matrix_free_apply_2d(wq_u_.B0, wq_v_.B1, x_ctrl, false);
        J21 = matrix_free_apply_2d(wq_u_.B1, wq_v_.B0, y_ctrl, false);
        J22 = matrix_free_apply_2d(wq_u_.B0, wq_v_.B1, y_ctrl, false);
    } else {
        Eigen::VectorXd wx = w_ctrl.cwiseProduct(x_ctrl);
        Eigen::VectorXd wy = w_ctrl.cwiseProduct(y_ctrl);

        Eigen::VectorXd W    = matrix_free_apply_2d(wq_u_.B0, wq_v_.B0, w_ctrl, false);
        Eigen::VectorXd dWdu = matrix_free_apply_2d(wq_u_.B1, wq_v_.B0, w_ctrl, false);
        Eigen::VectorXd dWdv = matrix_free_apply_2d(wq_u_.B0, wq_v_.B1, w_ctrl, false);

        Eigen::VectorXd Xnum    = matrix_free_apply_2d(wq_u_.B0, wq_v_.B0, wx, false);
        Eigen::VectorXd dXnumdu = matrix_free_apply_2d(wq_u_.B1, wq_v_.B0, wx, false);
        Eigen::VectorXd dXnumdv = matrix_free_apply_2d(wq_u_.B0, wq_v_.B1, wx, false);

        Eigen::VectorXd Ynum    = matrix_free_apply_2d(wq_u_.B0, wq_v_.B0, wy, false);
        Eigen::VectorXd dYnumdu = matrix_free_apply_2d(wq_u_.B1, wq_v_.B0, wy, false);
        Eigen::VectorXd dYnumdv = matrix_free_apply_2d(wq_u_.B0, wq_v_.B1, wy, false);

        Eigen::VectorXd Winv2 = W.cwiseInverse().cwiseAbs2();

        J11 = (dXnumdu.cwiseProduct(W) - Xnum.cwiseProduct(dWdu)).cwiseProduct(Winv2);
        J12 = (dXnumdv.cwiseProduct(W) - Xnum.cwiseProduct(dWdv)).cwiseProduct(Winv2);
        J21 = (dYnumdu.cwiseProduct(W) - Ynum.cwiseProduct(dWdu)).cwiseProduct(Winv2);
        J22 = (dYnumdv.cwiseProduct(W) - Ynum.cwiseProduct(dWdv)).cwiseProduct(Winv2);

        // Same quotient-rule ingredients, kept for the DISPLACEMENT field's
        // own rational correction in apply() (see the class docstring's
        // "NURBS RATIONALITY" note): w_proj_ = 1/W(xi), z_proj_ = grad(W)/W^2.
        w_ctrl_ = w_ctrl;
        w_proj_ = W.cwiseInverse();
        z_proj_u_ = dWdu.cwiseProduct(Winv2);
        z_proj_v_ = dWdv.cwiseProduct(Winv2);
    }
    rational_ = rational;

    const int nq = nq_u_ * nq_v_;
    Eigen::VectorXd detJ = J11.cwiseProduct(J22) - J12.cwiseProduct(J21);
    Eigen::VectorXd absDetJ = detJ.cwiseAbs();

    // invJ[alpha][l] = d(param_alpha)/d(phys_l), same convention as
    // PatchIntegrator's invJ11=J22/detJ etc.
    Eigen::VectorXd invJ[2][2];
    invJ[0][0] = J22.cwiseQuotient(detJ);
    invJ[0][1] = (-J12).cwiseQuotient(detJ);
    invJ[1][0] = (-J21).cwiseQuotient(detJ);
    invJ[1][1] = J11.cwiseQuotient(detJ);

    auto [lambda, mu] = extract_lame(law);

    // stiffness_property[I,J,alpha,beta] = |detJ| * [
    //     lambda*invJ[alpha][I]*invJ[beta][J]
    //   + mu*invJ[alpha][J]*invJ[beta][I]
    //   + mu*delta_IJ*(invJ[alpha][0]*invJ[beta][0] + invJ[alpha][1]*invJ[beta][1]) ]
    for (int I = 0; I < 2; ++I) {
        for (int J = 0; J < 2; ++J) {
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    Eigen::VectorXd term =
                        lambda * invJ[a][I].cwiseProduct(invJ[b][J])
                        + mu * invJ[a][J].cwiseProduct(invJ[b][I]);
                    if (I == J) {
                        term += mu * (invJ[a][0].cwiseProduct(invJ[b][0])
                                    + invJ[a][1].cwiseProduct(invJ[b][1]));
                    }
                    stiffness_property_[I][J][a][b] = term.cwiseProduct(absDetJ);
                }
            }
        }
    }
    (void)nq;
}

namespace {
const Eigen::SparseMatrix<double>& select_basis(const WeightedQuadrature1D& wq, int order) {
    return order == 0 ? wq.B0 : wq.B1;
}
const Eigen::SparseMatrix<double>& select_weight(const WeightedQuadrature1D& wq, int test_order, int trial_order) {
    int idx = 2 * test_order + trial_order;
    switch (idx) {
        case 0: return wq.W00;
        case 1: return wq.W01;
        case 2: return wq.W10;
        default: return wq.W11;
    }
}
} // namespace

Eigen::VectorXd WQMatrixFreeStiffness::apply(const Eigen::VectorXd& v) const {
    const int nctrl = nbctrlpts_u_ * nbctrlpts_v_;

    Eigen::VectorXd v_comp[2];
    for (int I = 0; I < 2; ++I) {
        v_comp[I].resize(nctrl);
        for (int k = 0; k < nctrl; ++k)
            v_comp[I](k) = v(patch_.dof_manager->get_global_dof(static_cast<size_t>(k), I));
    }

    Eigen::VectorXd out_comp[2] = {
        Eigen::VectorXd::Zero(nctrl), Eigen::VectorXd::Zero(nctrl)
    };

    // Value ("order 0 in both directions") gather/scatter, reused by the
    // NURBS correction terms below.
    const Eigen::SparseMatrix<double>& B0u = select_basis(wq_u_, 0);
    const Eigen::SparseMatrix<double>& B0v = select_basis(wq_v_, 0);
    const Eigen::SparseMatrix<double>& Wval_u = select_weight(wq_u_, 0, 0);
    const Eigen::SparseMatrix<double>& Wval_v = select_weight(wq_v_, 0, 0);
    Eigen::VectorXd w_proj_sq;
    if (rational_) w_proj_sq = w_proj_.cwiseProduct(w_proj_);
    auto z_proj = [this](int dir) -> const Eigen::VectorXd& {
        return dir == 0 ? z_proj_u_ : z_proj_v_;
    };

    for (int I = 0; I < 2; ++I) {
        for (int J = 0; J < 2; ++J) {
            // Trial-side field: raw DOFs for a B-spline patch, or the
            // NURBS-weighted DOFs (w_a * d_a) for a rational one -- see the
            // class docstring's "NURBS RATIONALITY" note. Interpolating this
            // "weighted field" with the raw B-spline basis and later dividing
            // by W(xi) (the w_proj_/z_proj_ terms below) reproduces the
            // correct rational basis function derivative via the quotient
            // rule, exactly as the constructor already does for the geometry
            // map.
            Eigen::VectorXd field_in = rational_
                ? Eigen::VectorXd(v_comp[J].cwiseProduct(w_ctrl_))
                : v_comp[J];

            for (int beta = 0; beta < 2; ++beta) {
                const Eigen::SparseMatrix<double>& Bu = select_basis(wq_u_, beta == 0 ? 1 : 0);
                const Eigen::SparseMatrix<double>& Bv = select_basis(wq_v_, beta == 1 ? 1 : 0);
                Eigen::VectorXd array_tmp = matrix_free_apply_2d(Bu, Bv, field_in, false);

                for (int alpha = 0; alpha < 2; ++alpha) {
                    Eigen::VectorXd coeff = stiffness_property_[I][J][alpha][beta];
                    if (rational_) coeff = coeff.cwiseProduct(w_proj_sq);
                    Eigen::VectorXd weighted = array_tmp.cwiseProduct(coeff);

                    const Eigen::SparseMatrix<double>& Wu =
                        select_weight(wq_u_, alpha == 0 ? 1 : 0, beta == 0 ? 1 : 0);
                    const Eigen::SparseMatrix<double>& Wv =
                        select_weight(wq_v_, alpha == 1 ? 1 : 0, beta == 1 ? 1 : 0);

                    out_comp[I] += matrix_free_apply_2d(Wu, Wv, weighted, false);
                }
            }

            if (!rational_) continue;

            // Quotient-rule cross terms, direct port of pymfiga's
            // NurbsOperations.compute_mf_scalar_gradu_gradv (see the class
            // docstring). The "clean gradient x gradient" term above only
            // captured the leading w_proj_^2 = 1/W^2 factor; these three
            // terms are what R_a = w_a*N_a/W(xi)'s own rationality (test
            // side) and the trial field's rationality contribute beyond
            // that, via z_proj_ = grad(W)/W^2.
            Eigen::VectorXd val_in = matrix_free_apply_2d(B0u, B0v, field_in, false);

            // Term 2 (value . value): coeff2 = sum_{a,b} SP[I,J,a,b]*z[a]*z[b]
            Eigen::VectorXd coeff2 = Eigen::VectorXd::Zero(val_in.size());
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    coeff2 += stiffness_property_[I][J][a][b].cwiseProduct(z_proj(a)).cwiseProduct(z_proj(b));
            out_comp[I] += matrix_free_apply_2d(Wval_u, Wval_v, val_in.cwiseProduct(coeff2), false);

            // Term 3 (subtract), test=value, trial=sum_beta grad(beta):
            // coeff3[beta] = sum_alpha SP[I,J,alpha,beta]*z[alpha]*w_proj
            Eigen::VectorXd term3_sum = Eigen::VectorXd::Zero(val_in.size());
            for (int beta = 0; beta < 2; ++beta) {
                Eigen::VectorXd coeff3 = Eigen::VectorXd::Zero(val_in.size());
                for (int a = 0; a < 2; ++a)
                    coeff3 += stiffness_property_[I][J][a][beta].cwiseProduct(z_proj(a)).cwiseProduct(w_proj_);
                const Eigen::SparseMatrix<double>& Bu = select_basis(wq_u_, beta == 0 ? 1 : 0);
                const Eigen::SparseMatrix<double>& Bv = select_basis(wq_v_, beta == 1 ? 1 : 0);
                term3_sum += matrix_free_apply_2d(Bu, Bv, field_in, false).cwiseProduct(coeff3);
            }
            out_comp[I] -= matrix_free_apply_2d(Wval_u, Wval_v, term3_sum, false);

            // Term 4 (subtract), test=grad(alpha), trial=value:
            // coeff4[alpha] = sum_beta SP[I,J,alpha,beta]*z[beta]*w_proj
            for (int alpha = 0; alpha < 2; ++alpha) {
                Eigen::VectorXd coeff4 = Eigen::VectorXd::Zero(val_in.size());
                for (int b = 0; b < 2; ++b)
                    coeff4 += stiffness_property_[I][J][alpha][b].cwiseProduct(z_proj(b)).cwiseProduct(w_proj_);
                const Eigen::SparseMatrix<double>& Wu = select_weight(wq_u_, alpha == 0 ? 1 : 0, 0);
                const Eigen::SparseMatrix<double>& Wv = select_weight(wq_v_, alpha == 1 ? 1 : 0, 0);
                out_comp[I] -= matrix_free_apply_2d(Wu, Wv, val_in.cwiseProduct(coeff4), false);
            }
        }
        // R_a = w_a*N_a/W(xi) -- the test side's own weight, applied once
        // after summing over J (it doesn't depend on J).
        if (rational_) out_comp[I] = out_comp[I].cwiseProduct(w_ctrl_);
    }

    Eigen::VectorXd out = Eigen::VectorXd::Zero(v.size());
    for (int I = 0; I < 2; ++I)
        for (int k = 0; k < nctrl; ++k)
            out(patch_.dof_manager->get_global_dof(static_cast<size_t>(k), I)) = out_comp[I](k);

    return out;
}
