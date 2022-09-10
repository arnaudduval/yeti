import numpy as np

qp_position = [0, 1, 4, 6, 3]
qp_position = np.array(qp_position)
qp_position.sort()
# qp_position = np.sort(qp_position)
print(qp_position)

p = [4, 5]

[a, b] = p
print(a)