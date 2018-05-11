import numpy as np
import polharmonic.shcoeffs as sh
import polharmonic.util as util
from scipy import special

class Illuminator:
    """An Illuminator is specified by its optical axis, numerical aperture, 
    the index of refraction of the sample, and polarizer orientation.

    By default we use the paraxial approximation.
    """
    def __init__(self, optical_axis=[0,0,1], na=0.8, n=1.33,
                 polarizer=None, paraxial=True, illuminate_all=False):
        self.optical_axis = optical_axis
        self.na = na
        self.n = n
        
        self.polarizer = polarizer
        self.paraxial = paraxial
        self.illuminate_all = illuminate_all

        self.alpha = np.arcsin(self.na/self.n)
        self.Atilde = (np.cos(self.alpha/2.0)**2)*np.cos(self.alpha)
        self.Btilde = (1.0/12.0)*(np.cos(self.alpha/2)**2 + 4*np.cos(self.alpha) + 7)

    def h(self):
        if self.illuminate_all:
            return sh.SHCoeffs([1.0, 0, 0, 0, 0, 0])
        
        p = self.polarizer
        if p is None:
            return sh.SHCoeffs([1, 0, 0, -self.Atilde/np.sqrt(5)])
        else:
            return sh.SHCoeffs([1,
                                -2*np.sqrt(0.6)*self.Btilde*p[0]*p[1],
                                0,
                                -self.Atilde/np.sqrt(5),
                                0,
                                np.sqrt(0.6)*self.Btilde*(p[0]**2 - p[1]**2)])

    # No spatioangular coupling -> h == H
    def H(self): 
        return self.h()
