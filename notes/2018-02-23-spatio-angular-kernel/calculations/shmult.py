import numpy as np
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import functools

# Analytical spherical fourier transform
import sympy.functions.special.spherical_harmonics as sh
def sh_mult(l1, m1, l2, m2):
    tex_string = ''
    for l in range(np.abs(l1 - l2), l1 + l2 + 2, 2):
        for m in range(-l, l+1):
            print("Integrating: "+ str(l) + ', ' + str(m))
            Znm1 = sh.Znm(l1, m1, theta, phi).expand(func=True)
            Znm2 = sh.Znm(l2, m2, theta, phi).expand(func=True)
            Znm = sh.Znm(l, m, theta, phi).expand(func=True)            
            theta_int = integrate(expand_trig(sin(theta)*Znm1*Znm2*Znm), (theta, 0, pi))
            final_int = integrate(expand_trig(theta_int.rewrite(cos)), (phi, 0, 2*pi))
            if final_int != 0:
                tex_string += latex(simplify(final_int)) + 'y_'+str(l)+'^{'+str(m)+'} + '
    print('y_'+str(l1)+'^{'+str(m1)+'}'+'y_'+str(l2)+'^{'+str(m2)+'} = ' + tex_string[:-2])            

# Polar fourier transform
print("Working...")

# Symbols
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)

# Calculate expansions
b = 1
for m1 in range(-b, b+1):
    for m2 in range(m1, b+1):
        sh_mult(b, m1, b, m2)

print('\n\n')
b = 2
for m1 in range(-b, b+1):
    for m2 in range(m1, b+1):
        sh_mult(b, m1, b, m2)
        
        
