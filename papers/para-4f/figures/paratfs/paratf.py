import numpy as np
np.seterr(all='ignore')
from scipy import special

def rect(r):
    return np.abs(r) < 0.5

def jinc(n, r):
    if n == 0:
        zero_val = np.pi/4
    else:
        zero_val = 0
    return np.where(r < 1e-5, zero_val, special.jn(n+1, np.pi*r)/(2*r))
        
def chat(n, nu):
    if n == 0:
        return np.where(nu < 1.0, 0.5*(np.arccos(nu) - nu*np.sqrt(1 - nu**2)), 0)
    elif n == 1:
        return np.where(nu < 1.0, 0.5*(np.arccos(nu) - nu*(3 - 2*(nu**2))*np.sqrt(1 - nu**2)), 0)
    else:
        return None

def Lambda(l):
    return np.sqrt(4*np.pi/(2*l + 1))

def norm(NA, n, nuc):
    return (12*nuc**2)/((np.pi*Lambda(0))*(2 + (NA/n)**2))
    
def dpsf(r, theta, NA=1.0, n=1.0, nuc=1.0):
    A = norm(NA, n, nuc)
    return A*(jinc(0, nuc*r)**2*np.sin(theta)**2 + (NA/n)**2*jinc(1, nuc*r)**2*np.cos(theta)**2)

def sdtf(nu, theta, NA=1.0, n=1.0, nuc=1.0):
    A = (norm(NA, n, nuc)/(nuc**2))
    return A*(chat(0, nu/nuc)*np.sin(theta)**2 + (NA/n)**2*chat(1, nu/nuc)*np.cos(theta)**2)

def adtf(r, l, NA=1.0, n=1.0, nuc=1.0):
    A = Lambda(l)*norm(NA, n, nuc)/3    
    if l == 0:
        return A*(2*jinc(0, nuc*r)**2 + (NA/n)**2*jinc(1, nuc*r)**2)
    elif l == 2:
        return A*(-2*jinc(0, nuc*r)**2 + 2*(NA/n)**2*jinc(1, nuc*r)**2)
    else:
        return 0

def satf(nu, l, NA=1.0, n=1.0, nuc=1.0):
    A = Lambda(l)*norm(NA, n, nuc)/3    
    if l == 0:
        return A*(2*chat(0, nu/nuc) + (NA/n)**2*chat(1, nu/nuc))
    elif l == 2:
        return A*(-2*chat(0, nu/nuc) + 2*(NA/n)**2*chat(1, nu/nuc))
    else:
        return 0
