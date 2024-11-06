from pysrc.lib.__init__ import *

# Create a 3-dimensional B-spline Curve
curve = NURBS.Curve()

# Set degree
curve.degree = 2

# Set control points (weights vector will be 1 by default)
# Use curve.ctrlptsw is if you are using homogeneous points as Pw
curve.ctrlpts = [[-1, 0, 0], [-1, 1, 0], [0, 1, 0]]

# Set knot vector
curve.knotvector = [0, 0, 0, 1, 1, 1]

# Set weights
curve.weights = [1, 1/np.sqrt(2), 1]

# Set evaluation delta (controls the number of curve points)
curve.delta = 0.05

# Get curve points (the curve will be automatically evaluated)
curve_points = np.array(curve.evalpts)

operations.insert_knot(curve, [0.5], [1])
print(curve.ctrlpts)
print(curve.knotvector)
print(curve.weights)

fig, ax = plt.subplots()
ax.plot(curve_points[:, 0], curve_points[:, 1])
ax.axis('equal')
fig.savefig("curveNURBS")