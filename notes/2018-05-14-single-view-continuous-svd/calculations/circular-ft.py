import numpy as np
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import functools

# Analytical spherical fourier transform

def cft(f, max_n=2):
    coeffs = []
    for n in range(-max_n, max_n+1):
        print("Integrating: "+ str(n))
        if n == 0:
            integrand = f/(2*sqrt(pi))
        if n > 0:
            integrand = f*cos(n*theta)/(sqrt(pi))
        if n < 0:
            integrand = f*sin(n*theta)/(sqrt(pi))
        theta_int = integrate(integrand, (theta, 0, 2*pi))            
        print(simplify(theta_int))
        coeffs.append(simplify(theta_int))
    return np.array(coeffs)

print("Working...")

# Symbols
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
a = Symbol('a', real=True)


I = 2*cos(phi - theta)**2 - 1
x = cft(I, max_n=2)

