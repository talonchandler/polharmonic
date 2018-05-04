# Tools for multiplying real spherical harmonics

# Based on:
# Herbert H.H. Homeier, E.Otto Steinborn,
# Some properties of the coupling coefficients of real spherical harmonics
# and their relation to Gaunt coefficients,
# Journal of Molecular Structure
# Volume 368, 1996, Pages 31-37, ISSN 0166-1280,
# https://doi.org/10.1016/S0166-1280(96)90531-X.

from util import *
import numpy as np
from sympy import *
from sympy.physics.wigner import gaunt, wigner_3j, clebsch_gordan
kd = KroneckerDelta

# Heaviside function
def hv(x):
    if x > 0:
        return 1
    else:
        return 0

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
def Rgaunt(l1, l2, l3, m1, m2, m3):
    result = 0
    for m1p in range(-l1, l1+1):
        U1 = U(l1, m1p, m1)
        for m2p in range(-l2, l2+1):
            U2 = U(l2, m2p, m2)
            for m3p in range(-l3, l3+1):
                U3 = U(l3, m3p, m3)
                result += U1*U2*U3*gaunt(l1, l2, l3, m1p, m2p, m3p)
    return sqrt(pi)*result

# Multiply two vectors of real spherical harmonic coefficients. Vectors must
# have the same length and be ordered like
#
# y_0^0, y_1^-1, y_1^0, y_1^1, y_2^-2, y_2^-1, ...
#
# Example:
#
# multiply_sh_coefficients([2, 1, 0, 0], [0, 1, 0, 0])
#
# gives the coefficients of (2y_0^0 + y_1^-1) x y_1^-1.
def multiply_sh_coefficients(a, b):
    maxl, m = i2lm(len(a))
    c = [0]*(maxl2maxi(2*(maxl-1)) + 1)
    for i, ai in enumerate(a):
        l1, m1 = i2lm(i)
        for j, bi in enumerate(b):
            l2, m2 = i2lm(j)
            for k, ci in enumerate(c):
                l3, m3 = i2lm(k)
                if ai != 0 and bi != 0:
                    c[k] += ai*bi*Rgaunt(l1, l2, l3, m1, m2, m3)

    return c

a = Symbol('a', real=True, positive=True)
b = Symbol('b', real=True, positive=True)
A = Symbol('A', real=True, positive=True)
B = Symbol('B', real=True, positive=True)
phi = Symbol('phi', real=True, positive=True)

# Polarized epi-illumination
# Unpolarized detection
result = multiply_sh_coefficients([a + 2*b, 0, 0, 0, 0, 0, (-a + 4*b)/sqrt(5), 0, 0],
                                  [1, 0, 0, 0, B*(sqrt(3)/sqrt(5))*(cos(phi)**2 - sin(phi)**2), 0, -A/sqrt(5), 0, -2*B*(sqrt(3)/sqrt(5))*cos(phi)*sin(phi)])

def rewrite(result):
    for res in result:
        print(collect(collect(expand(res), a), b))

import pdb; pdb.set_trace() 
print(simplify(result))
