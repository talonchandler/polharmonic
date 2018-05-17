from polharmonic import det, shcoeffs, gaunt
from sympy import *

At = Symbol('At')
a = Symbol('a')
b = Symbol('b')
z2 = Symbol('z2')
z_2 = Symbol('z_2')

# hlmn convention
# Symbolic
# h000e = Symbol('h000e')
# h200e = Symbol('h200e')

# h000d = Symbol('h000d')
# h002d = Symbol('h002d')
# h00_2d = Symbol('h00_2d')
# h200d = Symbol('h200d')
# h202d = Symbol('h202d')
# h20_2d = Symbol('h20_2d')
# h222d = Symbol('h222d')
# h2_2_2d = Symbol('h2_2_2d')

# Concrete
h000e = 1
h200e = -At/sqrt(5)

h000d = a + 2*b
h002d = 2*sqrt(pi)*b*z2
h00_2d = 2*sqrt(pi)*b*z_2
h200d = (-a + 4*b)/sqrt(5)
h202d = 4*sqrt(pi/5)*b*z2
h20_2d = 4*sqrt(pi/5)*b*z_2
h222d = sqrt(3/5)*a
h2_2_2d = sqrt(3/5)*a

excitation = [h000e, 0, 0, h200e, 0, 0]
# n = 
detection0 = [h000d, 0, 0, h200d, 0, 0]
detection2 = [h002d, 0, 0, h202d, 0, h222d]
detection_2 = [h00_2d, h2_2_2d, 0, h20_2d, 0, 0]

h0 = gaunt.multiply_sh_coefficients(excitation, detection0, evaluate=False)
h2 = gaunt.multiply_sh_coefficients(excitation, detection2, evaluate=False)
h_2 = gaunt.multiply_sh_coefficients(excitation, detection_2, evaluate=False)

hh0 = [simplify(x*(2*sqrt(pi))) for x in h0]
hh2 = [simplify(x*(2*sqrt(pi))) for x in h2]
hh_2 = [simplify(x*(2*sqrt(pi))) for x in h_2]
import pdb; pdb.set_trace() 
print(xx)
