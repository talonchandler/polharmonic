import numpy as np
import polharmonic.shcoeffs as sh
import polharmonic.tfcoeffs as tf
import polharmonic.util as util
from scipy import special

class Detector:
    """A Detector is specified by its optical axis, numerical aperture, 
    the index of refraction of the sample, and precence of a polarizer.

    By default we use the paraxial approximation.
    """
    def __init__(self, optical_axis=[0,0,1], na=0.8, n=1.33,
                 polarizer=False, paraxial=True, detect_all=False):
        self.optical_axis = optical_axis
        self.na = na
        self.n = n
        self.polarizer = polarizer
        self.paraxial = paraxial
        self.detect_all = detect_all

    def h(self, r, phi=0):
        if self.detect_all:
            return tf.TFCoeffs([[1.0, 0, 0, 0, 0, 0], 6*[0], 6*[0]])

        a = self.a
        b = self.b
        n0 = [a(r) + 2*b(r), 0, 0, (-a(r) + 4*b(r))/np.sqrt(5), 0, 0]
        if self.polarizer:
            n_2 = [2*b(r)*np.sin(2*phi), -np.sqrt(0.6)*a(r), 0, (4.0/np.sqrt(5))*b(r)*np.sin(2*phi), 0, 0]
            n2 = [2*b(r)*np.cos(2*phi), 0, 0, (4.0/np.sqrt(5))*b(r)*np.cos(2*phi), 0, np.sqrt(0.6)*a(r)]
            return tf.TFCoeffs([n0, n_2, n2])
        else:
            return tf.TFCoeffs([n0, 6*[0], 6*[0]])

    def H(self, nu, phi_nu=0):
        if self.detect_all:
            return tf.TFCoeffs([[1.0, 0, 0, 0, 0, 0], 6*[0], 6*[0]])

        A = self.A
        B = self.B
        C = self.C
        n0 = [A(nu) + 2*B(nu), 0, 0, (-A(nu) + 4*B(nu))/np.sqrt(5), 0, 0]
        if self.polarizer:
            n_2 = [2*C(nu)*np.sin(2*phi_nu), -np.sqrt(0.6)*A(nu), 0, (4.0/np.sqrt(5))*C(nu)*np.sin(2*phi_nu), 0, 0]
            n2 = [2*C(nu)*np.cos(2*phi_nu), 0, 0, (4.0/np.sqrt(5))*C(nu)*np.cos(2*phi_nu), 0, np.sqrt(0.6)*A(nu)]
            return tf.TFCoeffs([n0, n_2, n2])#/n0[0]
        else:
            return tf.TFCoeffs([n0, 6*[0], 6*[0]])#/n0[0]

        
    # PSF helper functions
    def a(self, r):
        if r == 0:
            return 1.0
        else:
            return (special.jn(1,2*np.pi*r)/(np.pi*r))**2

    def b(self, r):
        if r == 0:
            return 0
        else:
            return ((self.na/self.n)*(special.jn(2,2*np.pi*r)/(np.pi*r)))**2

    # OTF helper functions
    def myacos(self, r):
        if isinstance(r, np.ndarray):
            r[np.abs(r) < 2] = 2
        else:
            r = r if abs(r) < 2 else 2
        return np.arccos(np.abs(r/2))

    def mysqrt(self, r):
        if isinstance(r, np.ndarray):
            r[np.abs(r) < 2] = 2
        else:
            r = r if abs(r) < 2 else 2
        return (np.abs(r/2))*np.sqrt(1 - (np.abs(r/2))**2)

    def A(self, r):
        return (2/np.pi)*(self.myacos(r) - self.mysqrt(r))

    def B(self, r):
        N = (1.0/np.pi)*((self.na/self.n)**2)
        poly = (3.0 - 2.0*(np.abs(r/2)**2))
        return N*(self.myacos(r) - poly*self.mysqrt(r))

    def C(self, r):
        N = (1.0/np.pi)*((self.na/self.n)**2)
        poly  = (4.0/3.0)*(1.0 - 1.0*(np.abs(r/2)**2))
        return -N*poly*self.mysqrt(r)
