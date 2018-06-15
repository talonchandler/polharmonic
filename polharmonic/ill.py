import numpy as np
import polharmonic.shcoeffs as sh
import polharmonic.tfcoeffs as tf
import polharmonic.util as util
from scipy import special

class Illuminator:
    """An Illuminator is specified by its optical axis, numerical aperture, 
    the index of refraction of the sample, and polarizer orientation.

    By default we use the paraxial approximation.
    """
    def __init__(self, optical_axis=[0,0,1], na=0.8, n=1.33,
                 polarizer=True, paraxial=True, illuminate_all=False):
        self.optical_axis = optical_axis
        self.na = na
        self.n = n
        
        self.polarizer = polarizer
        self.paraxial = paraxial
        self.illuminate_all = illuminate_all

    def h(self):
        if self.illuminate_all:
            return tf.TFCoeffs([[1.0, 0, 0, 0, 0, 0], 6*[0], 6*[0]])

        n0 = [1 + (self.na/self.n)**2, 0, 0, -(1 - 2*((self.na/self.n)**2))/np.sqrt(5), 0, 0]
        n_2 = 6*[0]
        n2 = 6*[0]
        if self.polarizer:
            n_2 = [0, -np.sqrt(3/5), 0, 0, 0, 0]
            n2 = [0, 0, 0, 0, 0, np.sqrt(3/5)]
        if self.optical_axis == [1,0,0]: # x-illumination
            n0 = sh.SHCoeffs(n0).rotate().coeffs
            n_2 = sh.SHCoeffs(n_2).rotate().coeffs            
            n2 = sh.SHCoeffs(n2).rotate().coeffs

        return tf.TFCoeffs([n0, n_2, n2])

    # No detection spatioangular coupling -> h == H
    def H(self): 
        return self.h()
