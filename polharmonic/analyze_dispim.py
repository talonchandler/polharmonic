import numpy as np
from util import *
from dispim_coeffs import *
np.set_printoptions(precision=2)

# Globals
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
phi_pol = Symbol('phi_pol', real=True)
A = Symbol('A', real=True)
B = Symbol('B', real=True)

# Main script
print("Working...")

# Testing forward and inverse transforms
#prf = 1 - (sin(Theta)*cos(Phi)*sin(theta)*cos(phi) + sin(Theta)*sin(Phi)*sin(theta)*sin(phi) + cos(Theta)*cos(theta))**2
#prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(1 - (sin(theta)**2)*(cos(phi)**2)))
#prf = ((sin(theta)*sin(phi)*sin(phi_pol) - cos(theta)*cos(phi_pol))**2)*2*(A + B*(sin(theta)**2))
#prf = sin(theta)**2 + 2
#prf = sin(theta)**4*(cos(phi)**2)
# xxx = prf2Hlm(prf, max_l=4)
# yyy = Hlm2prf(xxx)
# print(simplify(prf - yyy.rewrite(cos)))

# Testing real forward and inverse transforms
#prf = 1 - (sin(Theta)*cos(Phi)*sin(theta)*cos(phi) + sin(Theta)*sin(Phi)*sin(theta)*sin(phi) + cos(Theta)*cos(theta))**2
#prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(1 - (sin(theta)**2)*(cos(phi)**2)))
#prf = ((sin(theta)*sin(phi)*sin(phi_pol) - cos(theta)*cos(phi_pol))**2)*2*(A + B*(sin(theta)**2))
#prf = sin(theta)**2 + 2
#prf = sin(theta)**4*(cos(phi)**2)
# xxx = prf2Hlm_real(prf, max_l=4)
# yyy = Hlm2prf_real(xxx)
# print(simplify(prf - yyy.rewrite(cos)))
# import pdb; pdb.set_trace() 

# Calculate forward model (time consuming)
# psi, labels = dispim_sys_matrix(0.8, 0.8, 1.33, coeff_file=None)
# dill.dump([psi, labels], open('psi.dat','wb'))

# Calculate forward model real (time consuming)
# psi, labels = dispim_sys_matrix_real(0.8, 0.8, 1.33)
# dill.dump([psi, labels], open('psi_dispim.dat','wb'))

# Load
f = open('psi_dispim.dat', 'rb')
psi, labels = dill.load(f)

M = psi.shape[0]
N = psi.shape[1]

ct = np.dot(psi.T, psi)
pinv = np.linalg.pinv(psi)
U, s, Vh = np.linalg.svd(psi.real)
import pdb; pdb.set_trace() 

# Check that the SVD is correct
# S = np.zeros((M, N), dtype=complex)
# S[:M, :M] = np.diag(s)
# print(np.allclose(psi, np.dot(U, np.dot(S, V))))
V = Vh.conj().T

# Plot
#plot_spherical([1.0], ['4,0'], filename='test.png')
for i in range(9):
    plot_spherical_real(V[:,i], labels, filename='spherical'+str(i)+'.png')
