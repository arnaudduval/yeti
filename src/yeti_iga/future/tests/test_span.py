import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "build"))


import numpy as np
from bspline import BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch
import matplotlib
import matplotlib.pyplot as plt

U = np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])
V = np.array([0., 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1., 1.])
W = np.array([0., 0., 0., 0., 0.3, 0.4, 0.8, 0.9, 1., 1., 1., 1.])

def test_BSpline_getters():
    b1 = BSpline(2, U)
    assert (b1.knot_vector == np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])).all()
    assert b1.degree == 2

def test_ND_BSpline():
    b1 = BSpline(2, U)
    b2 = BSpline(2, V)
    b3 = BSpline(3, W)
    surf = BSplineSurface(b1, b2)
    vol = BSplineVolume(b1, b2, b3)

    u = np.array([0.3, 0.5])
    assert (surf.find_span_nd(u) == [2, 4]).all()

    # TODO : ajouter un test de calcul des fonctions de base
    # print(surf.basis_funs_nd(span, u))

    u = np.array([0.3, 0.7, 0.5])
    assert(vol.find_span_nd(u) == [2, 5, 5]).all()

    # TODO : ajouter un test de calcul des fonctions de base
    # print(vol.basis_funs_nd(span, u).shape)

test_BSpline_getters()
test_ND_BSpline()




u = np.arange(0, 1, 0.1)








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



print("-------")

mgr = ControlPointManager(dim=2)

id0 = mgr.add_point([0.0, 0.0])
id1 = mgr.add_point([1.0, 0.0])
id2 = mgr.add_point([0.0, 1.0])
id3 = mgr.add_point([1.0, 1.0])
# TODO : les ajouter dans un autre ordre et avoir un mapping non ordonné


print("Nombre total de points :", mgr.n_points) # TODO : retourne un pointeur vers la fonction et non la valeur de retour
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
# ERROR : IL Y A UNE INVERSION DE DIRECTION DANS CE QUI EST FAIT CI-DESSOUS !!!
u = np.array([0.25, 0.75], dtype=np.float64)
spans = surf.find_span_nd(u)
print(spans)
tensor_basis = surf.basis_funs_nd(spans, u)
print(f'{tensor_basis = }')

local_pts = patch.local_control_point_view()[mapping, :]
dim_u, dim_v = tensor_basis.shape
local_pts_reshaped = local_pts.reshape((dim_u, dim_v, mgr.dim_phys))
surface_pt = np.tensordot(tensor_basis, local_pts_reshaped, axes=([0,1],[0,1]))

print(f"Surface evaluated at u={u[0]}, v={u[1]} : {surface_pt}")



print("======")
u = np.array([[0.25, 0.75],[0.7, 0.66]], dtype=np.float64)
spans = np.array([surf.find_span_nd(u[0, :]), surf.find_span_nd(u[1])])
print(spans)

print(patch.evaluate_patch_nd(spans, u))

