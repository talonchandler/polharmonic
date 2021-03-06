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
            coeffs.append(final_int)
    return np.array(coeffs)

# Polar fourier transform

print("Working...")

# Symbols
Theta = Symbol('Theta')
Phi = Symbol('Phi')
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)

# Calculate internal integrals

# Make these symbols for now
A = Symbol('A', real=True)
B = Symbol('B', real=True)
C = Symbol('C', real=True)
D = Symbol('D', real=True)

# Calculate intensity
I = A + B*(sin(theta)**2) + C*(sin(theta)**2)*cos(2*(phi - Phi))

# Expand onto spherical harmonics
x = sft(I, max_l=4, odd_l=False)
