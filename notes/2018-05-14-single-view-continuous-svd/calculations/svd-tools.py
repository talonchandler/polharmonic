from polharmonic import util as myutil
from sympy import *
from sympy.physics.wigner import gaunt, wigner_3j, clebsch_gordan
kd = KroneckerDelta
import numpy as np

# Illumination 
def hill(polarizer=True):
    At = Symbol('At')
    Bt = Symbol('Bt')
    
    n0 = [1, 0, 0, -At/sqrt(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if polarizer:
        n_2 = [0, -sqrt(3*pi/5)*Bt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        n2 = [0, 0, 0, 0, 0, sqrt(3*pi/5)*Bt, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        return Array([n0, n_2, n2])
    else:
        return Array([n0, 15*[0], 15*[0]])

# No detection spatioangular coupling -> h == H
def Hill(polarizer=True): 
    return hill(polarizer)

# Detection
def hdet(polarizer=False):
    a = Symbol('a')
    b = Symbol('b')
    phi = Symbol('phi')
    n0 = [a + 2*b, 0, 0, (-a + 4*b)/sqrt(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if polarizer:
        n_2 = [2*b*sin(2*phi), -sqrt(3/5)*a, 0, (4/sqrt(5))*b*sin(2*phi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        n2 = [2*b*cos(2*phi), 0, 0, (4.0/sqrt(5))*b*cos(2*phi), 0, sqrt(3/4)*a, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        return Array([n0, n_2, n2])
    else:
        return Array([n0, 15*[0], 15*[0]])

def Hdet(polarizer=False):
    A = Symbol('A')
    B = Symbol('B')
    C = Symbol('C')
    phi_nu = Symbol('phi_nu')
    n0 = [A + 2*B, 0, 0, (-A + 4*B)/sqrt(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if polarizer:
        n_2 = [2*C*cos(2*phi_nu), -sqrt(3/5)*A, 0, (4/sqrt(5))*C*cos(2*phi_nu), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        n2 = [2*C*sin(2*phi_nu), 0, 0, (4/sqrt(5))*C*sin(2*phi_nu), 0, sqrt(3/5)*A, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        return Array([n0, n_2, n2])
    else:
        return Array([n0, 15*[0], 15*[0]])

# Heaviside function
def hv(x):
    if x > 0:
        return 1
    else:
        return 0

# Calculate gaunt symbolic
# Unitary matrix that transforms complex sh to real sh
# See Eq. 12.
def U(l, m, mu):
    t1 = kd(m, 0)*kd(mu, 0)
    t2 = hv(mu)*kd(m, mu)
    t3 = hv(-mu)*I*((-1)**np.abs(m))*kd(m, mu)
    t4 = hv(-mu)*(-I)*kd(m, -mu)
    t5 = hv(mu)*((-1)**np.abs(m))*kd(m, -mu)
    return  t1 + ((t2 + t3 + t4 + t5)/sqrt(2))

# Real gaunt coefficients
# See Eqs. 26. The sympy gaunt function does not use a complex conjugate.
# This sum could be truncated using selection rules, but this is fairly quick.
def Rgaunt(l1, l2, l3, m1, m2, m3, evaluate=True):
    result = 0
    for m1p in range(-l1, l1+1):
        U1 = U(l1, m1p, m1)
        for m2p in range(-l2, l2+1):
            U2 = U(l2, m2p, m2)
            for m3p in range(-l3, l3+1):
                U3 = U(l3, m3p, m3)
                result += U1*U2*U3*gaunt(l1, l2, l3, m1p, m2p, m3p)
    if evaluate:
        return result.evalf()
    else:
        return result

# Compute and save an array with all of the gaunt coefficients up to specified
# band
def calc_gaunt_tensor(filename, lmax=4):
    jmax = myutil.maxl2maxj(lmax)
    G = np.zeros((jmax, jmax, jmax))
    GG = MutableDenseNDimArray(G)
    for index, g in np.ndenumerate(G):
        print(index)
        l1, m1 = myutil.j2lm(index[0])
        l2, m2 = myutil.j2lm(index[1])
        l3, m3 = myutil.j2lm(index[2])
        GG[index] = Rgaunt(l1, l2, l3, m1, m2, m3, evaluate=False)
    np.save(filename, GG)
    return GG

def multiply_tf_coeffs(in1, in2, P, G):
    out = np.zeros((P.shape[0], G.shape[0]))
    outA = MutableDenseNDimArray(out)
    for a in range(P.shape[0]):
        print(a)
        for b in range(P.shape[0]):
            for c in range(P.shape[0]):
                for d in range(G.shape[0]):
                    for e in range(G.shape[0]):
                        for f in range(G.shape[0]):
                            outA[c, f] += P[a, b, c]*G[d, e, f]*in1[a, d]*in2[b, e]
    return outA
                    
## MAIN ###

# Load/calculate tripling coefficients
import os
# G = calc_gaunt_tensor('gaunt_l4sym.npy', lmax=4)
G = Array(np.load(os.path.join(os.path.dirname(__file__), 'gaunt_l4sym.npy')), (15, 15, 15))

P = np.load(os.path.join(os.path.dirname(__file__), 'chcoeff_n2.npy'))
P = Array(P[:3, :3, :3])*sqrt(2*pi)
P = P.applyfunc(nsimplify)


# Calculate transfer function pieces
Hill_nopol = Hill(polarizer=False)
Hdet_pol = Hdet(polarizer=True)

Hill_pol = Hill(polarizer=True)
Hdet_nopol = Hdet(polarizer=False)

# Calculate complete transfer functions
H_poldet = multiply_tf_coeffs(Hill_nopol, Hdet_pol, P, G)
#H_polill = multiply_tf_coeffs(Hill_pol, Hdet_nopol, P, G)

# Calculate K
def calcK(H):
    K = MutableDenseNDimArray(np.zeros((P.shape[0], P.shape[0])))
    for i in range(P.shape[0]):
        for j in range(P.shape[0]):
            for k in range(G.shape[0]):
                K[i,j] += H[i,k]*H[j,k]
    return K

K = calcK(H_poldet)
K = K.tomatrix()
import pdb; pdb.set_trace() 
K.eigenvals()

import pdb; pdb.set_trace() 
