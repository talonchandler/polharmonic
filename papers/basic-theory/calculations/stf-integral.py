from sympy import *

a = Symbol('a')
A1 = Symbol('A1')
A2 = Symbol('A2')
theta = Symbol('theta')
phi = Symbol('phi')

Y00 = Ynm(0, 0, theta, phi).expand(func=True)
Y20 = Ynm(2, 0, theta, phi).expand(func=True)
Y40 = Ynm(4, 0, theta, phi).expand(func=True)

# Calc expansions
# hh = sin(theta)**2 + ((a**2)/2)*cos(theta)**2
hh = (sin(theta)**2 + ((a**2)/2)*cos(theta)**2)*(A1*sin(theta)**2 + A2*((a**2)/2)*cos(theta)**2)
c = 4*15*7/(2*sqrt(pi))
H0 = c*simplify(integrate(expand(2*pi*sin(theta)*hh*Y00), (theta, 0, pi)))
H2 = c*simplify(integrate(expand(2*pi*sin(theta)*hh*Y20), (theta, 0, pi)))
H4 = c*simplify(integrate(expand(2*pi*sin(theta)*hh*Y40), (theta, 0, pi)))
import pdb; pdb.set_trace() 
