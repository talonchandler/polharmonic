from sympy import *

# A = Symbol('A', real=True)
# B = Symbol('B', real=True)
alpha = Symbol('alpha', real=True)
n = Symbol('n', real=True)

A = (cos(alpha/2)**2)*(cos(alpha))
At = A.subs([(alpha, asin(n))])
As = simplify(expand(At))

B =  (cos(alpha)**2 + 4*cos(alpha) + 7)/12
Bt = B.subs([(alpha, asin(n))])
Bs = simplify(expand(Bt))
import pdb; pdb.set_trace() 
