""" In this file we test how to separate multivariate function as 
	product of univariate functions. Ex F(x, y) = G(x) H(y)
"""
import numpy as np

maxit = 10
Mex = np.array([[1, 4, 36], [16, 100, 1600], [9, 25, 1]])

n1, n2 = np.shape(Mex)
D1 = np.ones(n1); D2 = np.ones(n2)

for it in range(maxit):
	for i in range(n1):
		Div = Mex[i, :]/D2
		D1[i] = np.sqrt(np.max(Div)*np.min(Div))

	for j in range(n2):
		Div = Mex[:, j]/D1
		D2[j] = np.sqrt(np.max(Div)*np.min(Div))
		
Mapp = np.zeros(np.shape(Mex))
for i in range(n1):
    for j in range(n2):
        Mapp[i, j] = D1[i]*D2[j]

print(D1)
print(D2)
print(Mex)
print(Mapp)