import numpy as np
import dill
from util import *
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import sympy.functions.special.spherical_harmonics as sh

# Global initialize symbols
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
Theta = Symbol('Theta', real=True)
Phi = Symbol('Phi', real=True)

# Helper functions
def NAn2AB(NA, n):
    alpha = np.arcsin(NA/n)
    return alpha2A(alpha), alpha2B(alpha)

def alpha2A(alpha):
    ca = np.cos(alpha)
    return 0.25 - 0.375*ca + 0.125*(ca**3)

def alpha2B(alpha):
    ca = np.cos(alpha)
    return 0.1875*ca + 0.1875*(ca**3)

# Convert point response function to spherical harmonic coefficients
# Spherical Fourier transform
def prf2Hlm(prf, max_l=4):
    coeffs = []
    for l in range(0, max_l+2, 2):
        for m in range(0, l+1):
            print("Integrating: "+ str(l) + ', ' + str(m))            
            Ynm = sh.Ynm(l, m, theta, phi).expand(func=True)
            theta_int = integrate(expand(sin(theta)*conjugate(Ynm)*prf), (theta, 0, pi)) # theta integral
            final_int = simplify(integrate(expand_trig(theta_int), (phi, 0, 2*pi))) # phi integral
            coeffs.append((final_int, str(l)+','+str(m)))
    dtypes = [('expr', 'O'), ('lm', 'U3')]
    return np.array(coeffs, dtype=dtypes)

def prf2Hlm_real(prf, max_l=4):
    coeffs = []
    for l in range(0, max_l+2, 2):
        for m in range(-l, l+1):
            print("Integrating: "+ str(l) + ', ' + str(m))            
            Znm = sh.Znm(l, m, theta, phi).expand(func=True)
            theta_int = integrate(expand(sin(theta)*Znm*prf), (theta, 0, pi)) # theta integral
            final_int = simplify(integrate(expand_trig(theta_int), (phi, 0, 2*pi))) # phi integral
            coeffs.append((final_int, l, m))
    dtypes = [('expr', 'O'), ('l', 'i4'), ('m', 'i4')]
    return np.array(coeffs, dtype=dtypes)

# Convert spherical harmonic coefficients to prf
# Spherical Fourier transform inverse
def Hlm2prf(Hlm):
    prf = 0
    for term in Hlm:
        l = int(term['l'][0])
        m = int(term['m'][2])
        if m == 0:
            prf += term['expr']*sh.Ynm(l, m, theta, phi).expand(func=True)
        if m != 0:
            prf += (1/np.sqrt(2))*term['expr']*sh.Ynm(l, m, theta, phi).expand(func=True)
            prf += (1/np.sqrt(2))*((-1)**m)*conjugate(term['expr'])*sh.Ynm(l, -m, theta, phi).expand(func=True)
    return prf

def Hlm2prf_real(Hlm):
    prf = 0
    for term in Hlm:
        prf += term['expr']*sh.Znm(term['l'], term['m'], theta, phi).expand(func=True)
    return prf

# Coefficients for several geometries
# z light-sheet illumination, x widefield detection
def illzdetx(phi_pol, A, B):
    prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(1 - (sin(theta)**2)*(cos(phi)**2)))
    return prf2Hlm(prf)

def illzdetx_real(phi_pol, A, B):
    prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(1 - (sin(theta)**2)*(cos(phi)**2)))
    return prf2Hlm_real(prf)

# x light-sheet illumination, z widefield detection
def illxdetz(phi_pol, A, B):
    prf = ((sin(theta)*sin(phi)*sin(phi_pol) - cos(theta)*cos(phi_pol))**2)*2*(A + B*(sin(theta)**2))
    return prf2Hlm(prf)

def illxdetz_real(phi_pol, A, B):
    prf = ((sin(theta)*sin(phi)*sin(phi_pol) - cos(theta)*cos(phi_pol))**2)*2*(A + B*(sin(theta)**2))
    return prf2Hlm_real(prf)

def illzdetz_real(phi_pol, A, B):
    prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*sin(theta)**2)
    return prf2Hlm_real(prf)

# Single detector 
def illalldetpixel(Theta, Phi):
    prf = 1 - (sin(Theta)*cos(Phi)*sin(theta)*cos(phi) + sin(Theta)*sin(Phi)*sin(theta)*sin(phi) + cos(Theta)*cos(theta))**2
    return prf2Hlm(prf)

# Dispim system matrix
def dispim_sys_matrix(NAx, NAz, n):
    sys_matrix = []
    phi_pol = Symbol('phi_pol', real=True)
    
    Ax, Bx = NAn2AB(NAx, n)
    Az, Bz = NAn2AB(NAz, n)

    ixdz = illxdetz(phi_pol, Az, Bz)
    izdx = illzdetx(phi_pol, Ax, Bx)

    # Populate system matrix
    phi_pols=[0, np.pi/4, np.pi/2, 3*np.pi/4]
    for i, pol in enumerate(phi_pols):
        for j, (Hlm, lm) in enumerate(izdx):
            sys_matrix.append((np.complex(Hlm.evalf(subs={phi_pol: pol})), i, j, lm))
        for j, (Hlm, lm) in enumerate(ixdz):            
            sys_matrix.append((np.complex(Hlm.evalf(subs={phi_pol: pol})), i+len(phi_pols), j, lm))

    dtypes = [('Hlm', 'c16'), ('i', 'i4'), ('j', 'i4'), ('lm', 'U3')]
    
    arr = np.array(sys_matrix, dtype=dtypes)
    arr, labels = listmatrix2arraymatrix(arr)
    return arr, labels

def dispim_sys_matrix_real(NAx, NAz, n):
    sys_matrix = []
    phi_pol = Symbol('phi_pol', real=True)
    
    Ax, Bx = NAn2AB(NAx, n)
    Az, Bz = NAn2AB(NAz, n)

    ixdz = illxdetz_real(phi_pol, Az, Bz)
    izdx = illzdetx_real(phi_pol, Ax, Bx)

    # Populate system matrix
    phi_pols=[0, np.pi/4, np.pi/2, 3*np.pi/4]
    for i, pol in enumerate(phi_pols):
        for j, (Hlm, l, m) in enumerate(izdx):
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i, j, l, m))
        for j, (Hlm, l, m) in enumerate(ixdz):            
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i+len(phi_pols), j, l, m))

    dtypes = [('Hlm', 'c16'), ('i', 'i4'), ('j', 'i4'), ('l', 'i4'), ('m', 'i4')]
    
    arr = np.array(sys_matrix, dtype=dtypes)
    arr, labels = listmatrix2arraymatrix_real(arr)
    return arr, labels

def epi_sys_matrix(NA, n):
    sys_matrix = []
    phi_pol = Symbol('phi_pol', real=True)
    
    A, B = NAn2AB(NA, n)

    izdz = illzdetz_real(phi_pol, A, B)

    # Populate system matrix
    phi_pols=[0, np.pi/4, np.pi/2, 3*np.pi/4]
    for i, pol in enumerate(phi_pols):
        for j, (Hlm, l, m) in enumerate(izdz):
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i, j, l, m))

    dtypes = [('Hlm', 'c16'), ('i', 'i4'), ('j', 'i4'), ('l', 'i4'), ('m', 'i4')]
    
    arr = np.array(sys_matrix, dtype=dtypes)
    arr, labels = listmatrix2arraymatrix_real(arr)
    return arr, labels
