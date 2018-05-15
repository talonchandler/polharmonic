from polharmonic import det, shcoeffs, gaunt
from sympy import *

At = Symbol('At')
a = Symbol('a')
b = Symbol('b')
z1 = Symbol('z1')
z2 = Symbol('z2')

h00e = Symbol('h00e')
h20e = Symbol('h20e')
h00d = Symbol('h00d')
h20d = Symbol('h20d')
h2_2d = Symbol('h2_2d')
h22d = Symbol('h22d')

x = gaunt.multiply_sh_coefficients([h00e,0,0,h20e,0,0],
                                   [h00d,h2_2d,0,h20d,0,h22d],
                                   evaluate=False)
xx = [simplify(z*2*sqrt(pi)) for z in x]
print(xx)
