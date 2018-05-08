import numpy as np
import polharmonic.util as util
from scipy import special

class Detector:
    """A Detector is specified by its optical axis, numerical aperture, 
    the index of refraction of the sample, and polarizer orientation.

    By default we use the paraxial approximation.
    """
    def __init__(self, optical_axis=[0,0,1], na=0.8, n=1.33,
                 polarizer=None, paraxial=True):
        self.optical_axis = optical_axis        
        self.na = na
        self.n = n
        self.polarizer = polarizer
        self.paraxial = paraxial

    def h(self, r, phi=0):
        # TOOD: Arbitrary optical axes; relax paraxial approximation
        if self.polarizer is None:
            return [a(r) + 2*b(r), 0, 0, (-a(r) + 4*b(r))/np.sqrt(5)]
        else:
            r_vec = [r*cos(phi), r*sin(phi), 0]
            return [a(r) + 4*b(r)*np.dot(r_vec, self.polarizer)**2,
                    np.sqrt(0.6)*a(r)*(self.polarizer[0]**2 - self.polarizer[1]**2),
                    0,
                    (-a(r) + 8*b(r)*np.dot(r_vec, self.polarizer)**2)/np.sqrt(5),
                    0,
                    -2*np.sqrt(0.6)*a(r)*(self.polarizer[0]**2 - self.polarizer[1]**2)]

    def H(self, nu, phi_nu=0):
        # TOOD: Arbitrary optical axes; relax paraxial approximation
        if self.polarizer is None:
            return [a(r) + 2*b(r), 0, 0, (-a(r) + 4*b(r))/np.sqrt(5)]
        else:
            nu_vec = [r*cos(phi_nu), r*sin(phi_nu), 0]
            return [a(r) + 4*b(r)*np.dot(r_vec, self.polarizer)**2,
                    np.sqrt(0.6)*a(r)*(self.polarizer[0]**2 - self.polarizer[1]**2),
                    0,
                    (-a(r) + 8*b(r)*np.dot(r_vec, self.polarizer)**2)/np.sqrt(5),
                    0,
                    -2*np.sqrt(0.6)*a(r)*(self.polarizer[0]**2 - self.polarizer[1]**2)]

        
    # PSF helper functions
    def a(self, r):
        if r == 0:
            return 1.0
        else:
            return (special.jn(1,2*np.pi*r)/(np.pi*r))**2

    def b(self, r):
        if x == 0:
            return 0
        else:
            return ((self.NA/self.n)*(special.jn(2,2*np.pi*x)/(np.pi*x)))**2

    # OTF helper functions
    def myacos(self, r):
        return np.nan_to_num(np.arccos(np.abs(r/2)))

    def mysqrt(self, r):
        return np.nan_to_num((np.abs(r/2))*np.sqrt(1 - (np.abs(r/2))**2))

    def mycubert(self, r):
        return np.nan_to_num((np.abs(r/2))*((1 - (np.abs(r/2))**2)**(1.5)))

    def A(self, r):
        return (2/np.pi)*(myacos(r) - mysqrt(r))

    def B(self, r):
        N = (1.0/np.pi)*((self.NA/self.n)**2)
        poly = (3.0 - 2.0*(np.abs(r/2)**2))
        return N*(myacos(r) - poly*mysqrt(r))

    def C(self, r):
        N = (1.0/np.pi)*((self.NA/self.n)**2)
        poly  = (4.0/3.0)*(1.0 - 1.0*(np.abs(r/2)**2))
        return -N*poly*mysqrt(r)
        
    # def h00(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     if phi_p==None:
    #         return a(r_o) + 2*b(r_o, NA, n)
    #     else:
    #         return a(r_o) + 4*b(r_o, NA, n)*(np.cos(phi - phi_p)**2)

    # def h20(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     if phi_p==None:
    #         return (1/np.sqrt(5))*(-a(r_o) + 4*b(r_o, NA, n))
    #     else:
    #         return (1/np.sqrt(5))*(-a(r_o) + 8*b(r_o, NA, n)*(np.cos(phi - phi_p)**2))

    # def h22(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     if phi_p==None:
    #         return np.zeros(r_o.shape)
    #     else:
    #         return np.sqrt(3.0/5.0)*a(r_o)*(np.cos(phi_p)**2 - np.sin(phi_p)**2)

    # def h2_2(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     if phi_p==None:
    #         return np.zeros(r_o.shape)
    #     else:
    #         return -2*np.sqrt(3.0/5.0)*a(r_o)*np.cos(phi_p)*np.sin(phi_p)


    # def H00(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     N = (1 + (NA/n)**2)    
    #     if phi_p==None:
    #         return (A(x) + 2*B(x, NA=NA, n=n))/N
    #     else:
    #         return (A(x) + 2*B(x, NA=NA, n=n) + 2*C(x, NA=NA, n=n)*(np.cos(2*(phi-phi_p))))/N

    # def H20(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     N = np.sqrt(5)*(1 + (NA/n)**2)
    #     if phi_p==None:
    #         return (-A(x) + 4*B(x, NA=NA, n=n))/N
    #     else:
    #         return (-A(x) + 4*B(x, NA=NA, n=n) + 4*C(x, NA=NA, n=n)*(np.cos(2*(phi-phi_p))))/N

    # def H22(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     if phi_p==None:
    #         return np.zeros(x.shape)
    #     else:
    #         return np.sqrt(3.0/5.0)*(A(x)*(np.cos(phi_p)**2 - np.sin(phi_p)**2))/(1 + (NA/n)**2)

    # def H2_2(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    #     if phi_p==None:
    #         return np.zeros(x.shape)
    #     else:
    #         return -2*np.sqrt(3.0/5.0)*(A(x)*np.cos(phi_p)*np.sin(phi_p))/(1 + (NA/n)**2)

# def NAn2AB(NA, n):
#     alpha = np.arcsin(NA/n)
#     return alpha2A(alpha), alpha2B(alpha)

# def alpha2A(alpha):
#     ca = np.cos(alpha)
#     return 0.25 - 0.375*ca + 0.125*(ca**3)

# def alpha2B(alpha):
#     ca = np.cos(alpha)
#     return 0.1875*ca + 0.1875*(ca**3)
