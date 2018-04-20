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
            import pdb; pdb.set_trace() 
            theta_int = integrate(expand_trig(sin(theta)*Znm*f), (theta, 0, pi))
            final_int = integrate(expand_trig(theta_int.rewrite(cos)), (phi, 0, 2*pi))
            print(simplify(final_int))
            coeffs.append(final_int)
    return np.array(coeffs)

# Polar fourier transform

print("Working...")

# Symbols
# Theta = Symbol('Theta')
# Phi = Symbol('Phi')
theta = Symbol('theta')
phi = Symbol('phi', real=True)
x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)


s = Matrix([cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)])
R = Matrix([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
sp = R*s
thetap = acos(sp[2])
phip = atan(sp[1]/sp[0])
import pdb; pdb.set_trace() 

# Calculate intensity
I0  = simplify(sh.Znm(1, 1, thetap, phip).expand(func=True))
# Expand onto spherical harmonics
x = sft(I0, max_l=2, odd_l=False)
