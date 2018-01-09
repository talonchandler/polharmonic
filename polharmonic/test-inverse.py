import numpy as np
from util import *
from dispim_coeffs import *

# Globals
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
phi_pol = Symbol('phi_pol', real=True)
A = Symbol('A', real=True)
B = Symbol('B', real=True)

# Main script
print("Working...")

# Test the inverse
prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(1 - (sin(theta)**2)*(cos(phi)**2)))
#prf = 1 - (sin(theta)**2)*(sin(phi)**2)
#prf = sin(theta)**2
Hlm = prf2Hlm(prf)
prf2 = Hlm2prf(Hlm)
print(simplify(expand(prf).rewrite(cos) - expand(prf2).rewrite(cos)))
