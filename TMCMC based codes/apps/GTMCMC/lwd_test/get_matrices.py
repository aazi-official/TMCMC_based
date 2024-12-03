import numpy as np 

spls = np.loadtxt('samples.dat')
spls_cov = np.cov(spls,rowvar=0)
spls_mean = np.mean(spls, axis=0)
A = spls_cov

A_det = np.linalg.det(A)
if A_det < 1e-6:
	for i in range(A.shape[0]):
		for j in range(A.shape[0]):
			if i != j:
				A[i][j] = 0
A_inv = np.linalg.inv(A)
A_sqrt = np.linalg.cholesky(A)
A_det = np.linalg.det(A)				
np.savetxt('spls_mean.dat', spls_mean)
np.savetxt('spls_cov_sqrt.dat', A_sqrt)
np.savetxt('spls_cov_inv.dat', [A_det])
f=open('spls_cov_inv.dat','a')
np.savetxt(f, A_inv)