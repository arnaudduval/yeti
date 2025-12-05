#include "BSpline.hpp"

BSpline::BSpline(int degree_, py::array_t<double> kvector_) : degree(degree_) {
        // cpoy knotvector
        auto kv = kvector_.unchecked<1>();
        kvector.resize(kv.shape(0));
        for (ssize_t i = 0; i < kv.shape(0); ++i) kvector[i] = kv(i);
    }

py::array_t<double> BSpline::kvView() const {
    const ssize_t n = kvector.size();
    const double* data = kvector.data();

    // create dummy capsule to handle lifetime
    auto capsule = py::capsule((void*)data, [](void*) {
        // do nothing, memroy belongs to std::vector
    });

    return py::array_t<double>(
        {n},
        {sizeof(double)},
        data,
        capsule
    );
}

int BSpline::FindSpan(double u) const {
    int n = kvector.size() - degree - 2;
    if (u >= kvector[n+1]) return n;
    if (u <= kvector[degree]) return degree;

    int low = degree;
    int high = n+1;
    int mid = (low + high) / 2;

    while ((u < kvector[mid]) || u >= kvector[mid+1]) {
        if (u < kvector[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

py::array_t<double> BSpline::BasisFuns(int span, double u) const {
    int p = degree;

    // Crée un array Python de taille p+1
    py::array_t<double> result({p + 1});

    // Récupération d’un pointeur direct → zéro overhead
    double* out_ptr = result.mutable_data();

    // Appel direct à ta version interne optimisée
    BasisFuns_raw(span, u, out_ptr);

    return result;
}

void BSpline::BasisFuns_raw(int span, double u, double* out) const {
    int p = degree;
    const auto& U = kvector;

    std::vector<double> left(p+1), right(p+1);

    out[0] = 1.0;

    for (int j = 1; j <= p; ++j) {
        left[j]  = u - U[span + 1 - j];
        right[j] = U[span + j] - u;

        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            double temp = out[r] / (right[r+1] + left[j-r]);
            out[r] = saved + right[r+1] * temp;
            saved = left[j-r] * temp;
        }
        out[j] = saved;
    }
}

// Compute derivatives of basis functions
// cnvention : ders[k][i] -> k = 0..n, i = 0..p
// indexation : buffer[k*(p+1) + i]
// Algo Piegl & Tiller
void BSpline::BasisFunsDerivatives(int span, double u, int n, double* ders_out) const {
    const int p = degree;
    const std::vector<double>& U = kvector;

    if (n > p) n = p;

    // Offsets
    const int stride = p + 1;

    std::vector<double> ndu((p+1)*(p+1));
    std::vector<double> left(p+1), right(p+1);

    ndu[0] = 1.0;

    for (int j = 1; j <= p; ++j) {
        left[j] = u - U[span + 1 - j];
        right[j] = U[span + j] - u;

        double saved = 0.;
        for (int r = 0; r < j; ++r) {
            ndu[j*stride + r] = right[r+1] + left[j-r];
            double temp = ndu[r*stride + (j-1)] / ndu[j*stride + r];

            ndu[r*stride + j] = saved + right[r+1] * temp;
            saved = left[j-r] * temp;
        }
        ndu[j*stride + j] = saved;
    }

    // Load basis values
    for (int j = 0; j <= p; ++j)
        ders_out[j] = ndu[j*stride + p];

    // Compute derivatives
    std::vector<double> a(2*(p+1));
    for (int r = 0; r <= p; ++r) {
        int s1=0;
        int s2=1;
        a[0] = 1.;

        // Derivatives order k = 1..d
        for (int k = 1; k <= n; ++k) {
            double d = 0.;
            int rk = r - k;
            int pk = p - k;
            int j1, j2;

            if (r >= k) {
                a[s2*stride] = a[s1*stride] / ndu[(pk+1)*stride + rk];
                d = a[s2*stride] * ndu[rk*stride + pk];
            }
            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;
            if (r-1 <= pk)
                j2 = k-1;
            else
                j2 = p-r;

            for (int j = j1; j <= j2; ++j) {
                a[s2*stride+j] = (a[s1*stride+j] - a[s1*stride+j-1])/ndu[(pk+1)*stride+rk+j];
                d += a[s2*stride+j]*ndu[(rk+j)*stride + pk];
            }

            if (r <= pk) {
                a[s2*stride+k] = -a[s1*stride+k-1]/ndu[(pk+1)*stride+r];
                d += a[s2*stride + k] * ndu[r*stride+pk];
            }
            ders_out[k*stride+r] = d;
            int j = s1; s1=s2; s2=j;

        }
    }

    int r = p;
    for (int k = 1; k <= n; ++k)
    {
        for (int j = 0; j <= p; ++j)
            ders_out[k*stride + j] *= r;
        r *= (p-k);
    }
}


double BSpline::OneBasisFun(double u, int i) {
    int m = kvector.size() - 1;
    if (((i == 0) && (u == kvector[0])) || ( (i == m - degree - 1) && (u == kvector[m])))
        return 1.0;

    if ((u < kvector[i]) || (u >= kvector[i+degree+1]))
        return 0.0;

    std::vector<double> N(degree+1);
    for (int j = 0; j <= degree; j++) {
        if ((u >= kvector[i + j]) && u < kvector[i+j+1])
            N[j] = 1.0;
        else
            N[j] = 0.0;
    }

    double saved = 0.0;
    for (int k=1; k <= degree; k++) {
        if (N[0] == 0.0)
            saved = 0.0;
        else
            saved = ((u - kvector[i])*N[0]) / (kvector[i+k] - kvector[i]);

        for (int j = 0; j<degree-k+1; j++) {
            double left = kvector[i+j+1];
            double right = kvector[i+j+k+1];
            if (N[j+1] == 0.0) {
                N[j] = saved;
                saved = 0.0;
            }
            else {
                double temp = N[j+1] / (right - left);
                N[j] = saved + (right-u)*temp;
                saved = (u-left)*temp;
            }
        }

    }
    return N[0];
}

// Interface for Python binding
py::array_t<double> bspline_basis_funs_derivatives(const BSpline& self, int span, double u, int d) {
    const int p = self.getDegree();
    if (d > p) d = p;

    const ssize_t rows = d + 1;
    const ssize_t cols = p + 1;

    // output array
    py::array_t<double> arr({rows, cols});
    double* ptr = arr.mutable_data();

    self.BasisFunsDerivatives(span, u, d, ptr);

    return arr;
}