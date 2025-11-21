import numpy as np
from bspline import BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch
import matplotlib
import matplotlib.pyplot as plt

u = np.arange(0, 1, 0.1)
U = np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])

V = np.array([0., 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1., 1.])
W = np.array([0., 0., 0., 0., 0.3, 0.4, 0.8, 0.9, 1., 1., 1., 1.])


b1 = BSpline(2, U)
print(f'{b1.knot_vector = }')
print(f'{b1.degree = }')

# knots = np.arange(0., 1., 0.05)
# n = b1.knot_vector.shape[0] - b1.degree - 1
# print(f'{n=}')
# funs = np.empty((knots.shape[0], n))

# for j, k in enumerate(knots):
#     for i in range(0, n):
#         funs[j, i] = b1.one_basis_fun(k, i)

# matplotlib.use("Qt5Agg")  # or Qt5Agg if installed
# plt.plot(knots, funs)
# plt.show()


# Creation d'un espace 2D
b2 = BSpline(2, V)
b3 = BSpline(3, W)
surf = BSplineSurface(b1, b2)
vol = BSplineVolume(b1, b2, b3)

u = np.array([0.3, 0.5])
span = surf.find_span_nd(u)
print(surf.basis_funs_nd(span, u))

u = np.array([0.3, 0.7, 0.5])
span = vol.find_span_nd(u)
print(vol.basis_funs_nd(span, u).shape)

print("-------")

mgr = ControlPointManager(dim=2)

id0 = mgr.add_point([0.0, 0.0])
id1 = mgr.add_point([1.0, 0.0])
id2 = mgr.add_point([0.0, 1.0])
id3 = mgr.add_point([1.0, 1.0])

print("Nombre total de points :", mgr.n_points)
print("Coordonnées globales (zero-copy) :\n", mgr.coords_view())

U = np.array([0., 0., 1., 1.])
V = np.array([0., 0., 1., 1.])
su = BSpline(1, U)
sv = BSpline(1, V)
surf = BSplineSurface(su, sv)

mapping = np.array([0, 1, 2, 3], dtype=np.int64)
local_shape = [2, 2]
patch = Patch(surf, mgr, mapping.tolist(), local_shape)

# ---------------- Accès aux points locaux ----------------
print("\nPoints locaux via local_cp_view :")
for i in range(len(mapping)):
    pt_local = patch.local_cp_view(i)  # numpy array 1D (dim_phys,)
    print(f"Local point {i}: {pt_local}")

# ---------------- Vue sur tous les points globaux ----------------
all_pts = patch.local_control_point_view()
print("\nTous les points globaux (shape={}):\n{}".format(all_pts.shape, all_pts))

# ---------------- Exemple d'utilisation avec slicing ----------------
# Récupérer uniquement les points du patch local
local_pts = all_pts[mapping, :]
print("\nPoints du patch local via slicing numpy :\n", local_pts)

print("========")
u = np.array([0.25, 0.75], dtype=np.float64)
spans = surf.find_span_nd(u)
tensor_basis = surf.basis_funs_nd(spans, u)

local_pts = patch.local_control_point_view()[mapping, :]
dim_u, dim_v = tensor_basis.shape
local_pts_reshaped = local_pts.reshape((dim_u, dim_v, mgr.dim_phys))
surface_pt = np.tensordot(tensor_basis, local_pts_reshaped, axes=([0,1],[0,1]))

print(f"Surface evaluated at u={u[0]}, v={u[1]} : {surface_pt}")
