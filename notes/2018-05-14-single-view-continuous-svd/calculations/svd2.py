from sympy import *

# hlmn convention
# Symbolic
# H000 = Symbol('H000')
# H002 = Symbol('H002')
# H00_2 = Symbol('H00_2')
# H200 = Symbol('H200')
# H202 = Symbol('H202')
# H20_2 = Symbol('H20_2')
# H222 = Symbol('H222')
# H400 = Symbol('H400')
# H402 = Symbol('H402')
# H40_2 = Symbol('H40_2')
# H422 = Symbol('H422')

A = Symbol('A')
B = Symbol('B')
C = Symbol('C')
H000 = A + B
H002 = C
H00_2 = C
H200 = A + B
H202 = C
H20_2 = C
H222 = A
H400 = A - B
H402 = C
H40_2 = C
H422 = A

K00 = H00_2**2 + H222**2 + H20_2**2 + H422**2 + H40_2**2
K11 = H000**2 + H200**2 + H400**2
K22 = H002**2 + H222**2 + H202**2 + H422**2 + H402**2
K10 = H00_2*H000 + H20_2*H200 + H40_2*H400
K02 = H00_2*H002 + H20_2*H202 + H40_2*H402
K12 = H000*H002 + H200*H202 + H400*H402

K = Matrix([[K00, K10, K02],
           [K10, K11, K12],
           [K02, K12, K22]])

import pdb; pdb.set_trace() 
