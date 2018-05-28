import numpy as np
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import functools

# Analytical spherical fourier transform
import sympy.functions.special.spherical_harmonics as sh
def sft(f, max_l=4, odd_l=False):
    coeffs = []
    for l in range(0 + odd_l, max_l+2, 2):
        for m in range(-l, l+1):
            print("Integrating: "+ str(l) + ', ' + str(m))
            Znm = sh.Znm(l, m, theta, phi).expand(func=True)
            Znm = simplify(Znm) # Try this
            # import pdb; pdb.set_trace() 
            theta_int = integrate(expand_trig(sin(theta)*Znm*f), (theta, 0, pi))
            final_int = integrate(expand_trig(theta_int.rewrite(cos)), (phi, 0, 2*pi))
            print(simplify(final_int))
            coeffs.append(simplify(final_int))
    return np.array(coeffs)

# Polar fourier transform

print("Working...")

# Symbols
Theta = Symbol('Theta')
Phi = Symbol('Phi')
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
a = Symbol('a', real=True)
n = Symbol('n', real=True)

# Calculate internal integrals

# Make these symbols for now
A = Symbol('A', real=True)
B = Symbol('B', real=True)
C = Symbol('C', real=True)
D = Symbol('D', real=True)

# Calculate intensity
I = D*(A + B*(sin(theta)**2) + C*(sin(theta)**2)*cos(2*(phi - Phi)))

# Expand onto spherical harmonics
x = sft(I, max_l=2, odd_l=False)

x = [y.subs(A, (1/4) - (3/8)*cos(a) + (1/8)*(cos(a)**3)) for y in x]
x = [y.subs(B, (3/16)*cos(a) - (3/16)*(cos(a)**3)) for y in x]
x = [y.subs(C, (7/32) - (3/32)*cos(a) - (3/32)*(cos(a)**2) - (1/32)*(cos(a)**3)) for y in x]
x = [y.subs(D, 4/(3*(1 - cos(a)))) for y in x]

x = [nsimplify(y) for y in x]

xx = Array(x)
xx = ((3/(2*sqrt(pi)))*xx).applyfunc(simplify)

xxx = [simplify(y.subs(a, asin(n))) for y in xx]
import pdb; pdb.set_trace() 

