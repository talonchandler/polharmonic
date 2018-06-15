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
    def __init__(self, optical_axis=[0,0,1], na=0.8, n=1.33, sigma_ax=1.0,
                 polarizer=False, paraxial=True, detect_all=False):
        self.optical_axis = optical_axis
        self.na = na
        self.n = n
        self.polarizer = polarizer
        self.paraxial = paraxial
        self.detect_all = detect_all
        self.sigma_ax = sigma_ax

    def h(self, x, y, z):
        if self.detect_all:
            return tf.TFCoeffs([[1.0, 0, 0, 0, 0, 0], 6*[0], 6*[0]])

        # Find cylindrical coordinates based on view
        if self.optical_axis == [0,0,1]: # z-detection
            r = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            z_ax = z
        elif self.optical_axis == [1,0,0]: # x-detection
            r = np.sqrt(y**2 + z**2)
            phi = np.arctan2(z,y)
            z_ax = x

        a = self.a
        b = self.b
        n0 = [a(r) + 2*b(r), 0, 0, (-a(r) + 4*b(r))/np.sqrt(5), 0, 0]
        n_2 = 6*[0]
        n2 = 6*[0]
        if self.polarizer:
            n_2 = [2*b(r)*np.sin(2*phi), -np.sqrt(0.6)*a(r), 0, (4.0/np.sqrt(5))*b(r)*np.sin(2*phi), 0, 0]
            n2 = [2*b(r)*np.cos(2*phi), 0, 0, (4.0/np.sqrt(5))*b(r)*np.cos(2*phi), 0, np.sqrt(0.6)*a(r)]
        if self.optical_axis == [1,0,0]: # x-detection            
            n0 = sh.SHCoeffs(n0).rotate().coeffs
            n_2 = sh.SHCoeffs(n_2).rotate().coeffs            
            n2 = sh.SHCoeffs(n2).rotate().coeffs
            
        return tf.TFCoeffs([n0, n_2, n2])
    
    def H(self, x, y, z):
        if self.detect_all:
            return tf.TFCoeffs([[1.0, 0, 0, 0, 0, 0], 6*[0], 6*[0]])

        # Find cylindrical coordinates based on view
        if self.optical_axis == [0,0,1]: # z-detection
            nu = np.sqrt(x**2 + y**2)
            phi_nu = np.arctan2(y, x)
            z_ax = z
        elif self.optical_axis == [1,0,0]: # x-detection
            nu = np.sqrt(y**2 + z**2)
            phi_nu = np.arctan2(z, y)
            z_ax = x

        A = self.A
        B = self.B
        C = self.C
        n0 = [A(nu) + 2*B(nu), 0, 0, (-A(nu) + 4*B(nu))/np.sqrt(5), 0, 0]
        n_2 = 6*[0]
        n2 = 6*[0]
        if self.polarizer:
            n_2 = [2*C(nu)*np.sin(2*phi_nu), -np.sqrt(0.6)*A(nu), 0, (4.0/np.sqrt(5))*C(nu)*np.sin(2*phi_nu), 0, 0]
            n2 = [2*C(nu)*np.cos(2*phi_nu), 0, 0, (4.0/np.sqrt(5))*C(nu)*np.cos(2*phi_nu), 0, np.sqrt(0.6)*A(nu)]
        if self.optical_axis == [1,0,0]: # x-detection            
            n0 = sh.SHCoeffs(n0).rotate().coeffs
            n_2 = sh.SHCoeffs(n_2).rotate().coeffs            
            n2 = sh.SHCoeffs(n2).rotate().coeffs

        ax_weight = np.exp(-(z_ax**2)/(2*(self.sigma_ax**2)))/(np.sqrt(2*np.pi*(self.sigma_ax**2)))
        return tf.TFCoeffs([n0, n_2, n2])*ax_weight

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
