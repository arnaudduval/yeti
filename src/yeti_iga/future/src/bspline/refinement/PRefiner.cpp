#include <algorithm>
#include <stdexcept>
#include <pybind11/numpy.h>
#include "refinement/PRefiner.hpp"
#include "refinement/HRefiner.hpp"
#include "BSplineTensor.hpp"


// ---------------------------------------------------------------------------
// 1-D elevation (Piegl & Tiller A5.9) with coefficient tracking
// ---------------------------------------------------------------------------

void PRefiner::elevate_1d(int p,
                           const std::vector<double>& U,
                           Eigen::MatrixXd& T1d,
                           std::vector<double>& new_kv)
{
    int n   = static_cast<int>(U.size()) - p - 2;  // last old CP index
    int m   = n + p + 1;                            // last knot index
    int ph  = p + 1;                                // new degree
    int npts = n + 1;                               // number of original CPs

    // Bezier degree-elevation coefficients for t=1:
    //   bezalfs[i][i-1] = i/(p+1),   bezalfs[i][i] = (p+1-i)/(p+1)
    //   bezalfs[0][0] = bezalfs[ph][p] = 1
    std::vector<std::vector<double>> bezalfs(ph + 1, std::vector<double>(p + 1, 0.0));
    bezalfs[0][0]   = 1.0;
    bezalfs[ph][p]  = 1.0;
    for (int i = 1; i <= p; i++) {
        bezalfs[i][i-1] = static_cast<double>(i)   / ph;
        bezalfs[i][i]   = static_cast<double>(ph-i) / ph;
    }

    // Each "control point" is tracked as a coefficient row-vector of length npts.
    auto zero_row = [&]() { return Eigen::RowVectorXd::Zero(npts); };

    std::vector<Eigen::RowVectorXd> bpts    (p  + 1, zero_row());
    std::vector<Eigen::RowVectorXd> ebpts   (ph + 2, zero_row());
    std::vector<Eigen::RowVectorXd> Nextbpts(p  + 1, zero_row());
    std::vector<double> alphas(p + 1, 0.0);

    // Upper bounds for output sizes
    int max_knots = 2 * m + 4;
    std::vector<double>              Uh(max_knots + 1, 0.0);
    std::vector<Eigen::RowVectorXd>  Q (npts + m + 2, zero_row());

    // --- Initialise ---
    int mh   = ph;
    int kind = ph + 1;
    int r    = -1;
    int a    = p;
    int b    = p + 1;
    int cind = 1;
    double ua = U[0];

    Q[0](0) = 1.0;                                // Q[0] = P[0]
    for (int i = 0; i <= ph; i++) Uh[i] = ua;    // ph+1 leading knots

    for (int i = 0; i <= p; i++) bpts[i](i) = 1.0; // bpts[i] = P[i]

    // --- Main loop over knot spans ---
    while (b < m) {
        int i = b;
        while (b < m && U[b] == U[b+1]) b++;
        int    mult = b - i + 1;
        mh         = mh + mult + 1;
        double ub  = U[b];
        int    oldr = r;
        r = p - mult;

        int lbz = (oldr > 0) ? (oldr + 2) / 2 : 1;
        int rbz = (r > 0)    ? ph - (r + 1) / 2 : ph;

        // Knot insertion to prepare next Bezier segment
        if (r > 0) {
            double numer = ub - ua;
            for (int k = p; k >= mult + 1; k--)
                alphas[k - mult - 1] = numer / (U[a + k] - ua);
            for (int j = 1; j <= r; j++) {
                int save = r - j;
                int s    = mult + j;
                for (int k = p; k >= s; k--) {
                    double al = alphas[k - s];
                    bpts[k] = al * bpts[k] + (1.0 - al) * bpts[k-1];
                }
                Nextbpts[save] = bpts[p];
            }
        }

        // Elevate current Bezier segment
        for (int ii = lbz; ii <= ph; ii++) {
            ebpts[ii] = zero_row();
            int mpi = std::min(p, ii);
            for (int j = std::max(0, ii - 1); j <= mpi; j++)
                ebpts[ii] += bezalfs[ii][j] * bpts[j];
        }

        // Remove knot ua oldr times (only if oldr > 1)
        if (oldr > 1) {
            int    first = kind - 2;
            int    last  = kind;
            double den   = ub - ua;
            double bet   = (ub - Uh[kind-1]) / den;

            for (int tr = 1; tr <= oldr - 1; tr++) {
                int ii = first, jj = last;
                int kj = jj - kind + 1;
                while (jj - ii > tr) {
                    if (ii < cind) {
                        double alf = (ub - Uh[ii]) / (ua - Uh[ii]);
                        Q[ii] = alf * Q[ii] + (1.0 - alf) * Q[ii-1];
                    }
                    if (jj >= lbz) {
                        double coeff;
                        if (jj - tr <= kind - ph + oldr) {
                            coeff = (ub - Uh[jj - tr]) / den;
                        } else {
                            coeff = bet;
                        }
                        ebpts[kj] = coeff * ebpts[kj] + (1.0 - coeff) * ebpts[kj+1];
                    }
                    ii++; jj--; kj--;
                }
                first--; last++;
            }
        }

        // Store interior knots (skip for first span where a == p)
        if (a != p) {
            for (int ii = 0; ii < ph - oldr; ii++) {
                Uh[kind++] = ua;
            }
        }

        // Store elevated Bezier CPs into Q
        for (int j = lbz; j <= rbz; j++) {
            Q[cind++] = ebpts[j];
        }

        if (b < m) {
            // Slide window: blend saved CPs + next segment
            for (int j = 0; j < r; j++)              bpts[j] = Nextbpts[j];
            for (int j = r; j <= p; j++) {
                bpts[j] = zero_row();
                bpts[j](b - p + j) = 1.0;            // bpts[j] = P[b-p+j]
            }
            a = b; b++; ua = ub;
        } else {
            for (int ii = 0; ii <= ph; ii++) Uh[kind + ii] = ub;
        }
    }

    int nh = mh - ph - 1;   // last new CP index

    // Build outputs
    T1d.resize(nh + 1, npts);
    for (int i = 0; i <= nh; i++) T1d.row(i) = Q[i];

    new_kv.assign(Uh.begin(), Uh.begin() + mh + 1);
}


// ---------------------------------------------------------------------------
// PRefiner::refine  (tensor-product: apply elevate_1d along direction_)
// ---------------------------------------------------------------------------

void PRefiner::refine(Patch& patch, Eigen::MatrixXd& transition_matrix) const
{
    if (direction_ < 0 || direction_ >= static_cast<int>(patch.tensor.components.size()))
        throw std::invalid_argument("Invalid direction for degree elevation.");
    if (n_elevations_ <= 0)
        throw std::invalid_argument("n_elevations must be >= 1.");

    // Initialise composed T as identity, then fold n_elevations_ single-step elevations.
    size_t nb_cp_init = patch.global_indices.size();
    transition_matrix = Eigen::MatrixXd::Identity(nb_cp_init, nb_cp_init);

    for (int elev = 0; elev < n_elevations_; ++elev) {
        Eigen::MatrixXd T_step;
        apply_one_elevation(patch, T_step);
        transition_matrix = T_step * transition_matrix;
    }
}

// Single elevation step (factored out so the loop stays clean).
void PRefiner::apply_one_elevation(Patch& patch, Eigen::MatrixXd& transition_matrix) const
{
    const BSpline& spline = patch.tensor.components[direction_];
    int p = spline.getDegree();

    Eigen::MatrixXd T1d;
    std::vector<double> new_kv;
    elevate_1d(p, spline.getKnotVector(), T1d, new_kv);

    int n_old = static_cast<int>(patch.local_shape[direction_]);
    int n_new = static_cast<int>(T1d.rows());

    size_t ndim       = patch.tensor.components.size();
    size_t nb_old_cp  = patch.global_indices.size();
    size_t dim_phys   = patch.cp_manager->dim_phys;

    // New shape
    std::vector<size_t> new_local_shape(patch.local_shape);
    new_local_shape[direction_] = n_new;

    size_t nb_new_cp = 1;
    for (size_t d = 0; d < ndim; d++) nb_new_cp *= new_local_shape[d];

    // u-fastest strides: direction 0 fastest (stride[0]=1)
    std::vector<size_t> old_stride(ndim), new_stride(ndim);
    old_stride[0] = new_stride[0] = 1;
    for (size_t d = 1; d < ndim; d++) {
        old_stride[d] = old_stride[d-1] * patch.local_shape[d-1];
        new_stride[d] = new_stride[d-1] * new_local_shape[d-1];
    }

    size_t nb_lines = nb_old_cp / static_cast<size_t>(n_old);

    // Build transition matrix and compute new CP coords in one pass
    transition_matrix = Eigen::MatrixXd::Zero(nb_new_cp, nb_old_cp);
    std::vector<std::vector<double>> new_coords(nb_new_cp,
                                                std::vector<double>(dim_phys, 0.0));

    std::vector<size_t> other_idx(ndim, 0);

    for (size_t line = 0; line < nb_lines; line++) {
        size_t ls_old = 0, ls_new = 0;
        for (size_t d = 0; d < ndim; d++) {
            if (d != static_cast<size_t>(direction_)) {
                ls_old += other_idx[d] * old_stride[d];
                ls_new += other_idx[d] * new_stride[d];
            }
        }

        for (int ni = 0; ni < n_new; ni++) {
            size_t nf = ls_new + static_cast<size_t>(ni) * new_stride[direction_];
            for (int oi = 0; oi < n_old; oi++) {
                double coeff = T1d(ni, oi);
                if (std::abs(coeff) < 1e-15) continue;
                size_t of = ls_old + static_cast<size_t>(oi) * old_stride[direction_];
                transition_matrix(nf, of) = coeff;
                const double* src = patch.local_cp_ptr(of);
                for (size_t d = 0; d < dim_phys; d++)
                    new_coords[nf][d] += coeff * src[d];
            }
        }

        for (int d = static_cast<int>(ndim) - 1; d >= 0; d--) {
            if (d == direction_) continue;
            if (++other_idx[d] < patch.local_shape[d]) break;
            other_idx[d] = 0;
        }
    }

    // Replace cp_manager contents with the new CPs (sequential global indices)
    {
        std::lock_guard<std::mutex> lock(patch.cp_manager->mtx);
        patch.cp_manager->coords.resize(nb_new_cp * dim_phys);
        for (size_t i = 0; i < nb_new_cp; i++)
            for (size_t d = 0; d < dim_phys; d++)
                patch.cp_manager->coords[i * dim_phys + d] = new_coords[i][d];
    }

    // Update patch
    std::vector<size_t> new_global_indices(nb_new_cp);
    std::iota(new_global_indices.begin(), new_global_indices.end(), 0);
    patch.global_indices = std::move(new_global_indices);
    patch.local_shape    = std::move(new_local_shape);

    std::vector<BSpline> new_components = patch.tensor.components;
    py::array_t<double> py_new_kv = py::cast(new_kv);
    new_components[direction_] = BSpline(p + 1, py_new_kv);
    patch.tensor = BSplineTensor(new_components);

    if (patch.dof_manager) {
        patch.dof_manager = std::make_shared<PatchDOFManager>(
            *patch.dof_manager, patch.global_indices);
    }
}


// ---------------------------------------------------------------------------
// PRefiner::refine_1d  — fast path: no full nD matrix
// ---------------------------------------------------------------------------

void PRefiner::refine_1d(Patch& patch, Eigen::MatrixXd& T_1d) const
{
    if (direction_ < 0 || direction_ >= static_cast<int>(patch.tensor.components.size()))
        throw std::invalid_argument("Invalid direction for degree elevation.");

    // Track the evolving 1D spline locally — no CP update per step.
    BSpline spline_1d = patch.tensor.components[direction_];
    int n = static_cast<int>(patch.local_shape[direction_]);
    T_1d = Eigen::MatrixXd::Identity(n, n);

    for (int elev = 0; elev < n_elevations_; ++elev) {
        int p = spline_1d.getDegree();
        Eigen::MatrixXd T_step;
        std::vector<double> new_kv;
        elevate_1d(p, spline_1d.getKnotVector(), T_step, new_kv);
        spline_1d = BSpline(p + 1, py::cast(new_kv));
        T_1d = T_step * T_1d;
    }

    // Apply composed T_1d to patch CPs once.
    HRefiner::apply_1d_cp_update(patch, direction_, T_1d);

    std::vector<BSpline> new_components = patch.tensor.components;
    new_components[direction_] = spline_1d;
    patch.tensor = BSplineTensor(new_components);
}
