import numpy as np
from util import *
from dispim_coeffs import *
import scipy
np.set_printoptions(precision=2)

# Globals
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
phi_pol = Symbol('phi_pol', real=True)
A = Symbol('A', real=True)
B = Symbol('B', real=True)

# Main script
print("Working...")

# Test forward and inverse transforms
# prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*sin(theta)**2)
# Hlm = prf2Hlm_real(prf, max_l=4)
# prf2 = Hlm2prf_real(Hlm)
# print(simplify(prf - prf2.rewrite(cos)))
# import pdb; pdb.set_trace() 

# Calculate forward model (time consuming)
psi, labels = epi_sys_matrix(0.8, 1.33)
dill.dump([psi, labels], open('psi_epi.dat','wb'))

# Load
f = open('psi_epi.dat', 'rb')
psi, labels = dill.load(f)

M = psi.shape[0]
N = psi.shape[1]

ct = np.dot(psi.T, psi)
pinv = np.linalg.pinv(psi)
U, s, Vh = np.linalg.svd(psi.real)
V = Vh.conj().T

# Plot singular distributions
# plot_spherical_real([1.0], ['3,0'], filename='test.png')
import pdb; pdb.set_trace()
for i in range(9):
    plot_spherical_real(V[:,i], labels, filename='spherical'+str(i)+'.png')
