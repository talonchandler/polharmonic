from sympy import *

# l,m,n
At = Symbol('At')
A = Symbol('A')
Bt = Symbol('Bt')
B = Symbol('B')
H000 = (At/10 + 1/2)*A + (-2*At/5 + 1)*B
H222 = (3*sqrt(15)*Bt/35)*(3*A/2 + B)
H200 = (-sqrt(5)*At/14 + sqrt(5)/10)*A + (-11*sqrt(5)*At/35 + 2/sqrt(5))*B
H422 = 6*sqrt(5)*Bt/35*(-A/4 + B)
H400 = 3*At/35*(A - 4*B)

m0 = H000**2 + H200**2 + H400**2
m1 = H222**2 + H422**2

import pdb; pdb.set_trace() 
