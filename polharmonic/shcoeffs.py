import numpy as np
from polharmonic import util, gaunt
import matplotlib.pyplot as plt

class SHCoeffs:
    """An SHCoeffs object stores real spherical harmonic coefficients for even
    bands. It provides methods for adding and multiplying these coefficients.

    Inspired by a similar class in SHTOOLS https://shtools.github.io/SHTOOLS/

    """

    def __init__(self, coeffs):
        self.lmax, mm = util.j2lm(len(coeffs) - 1)
        self.jmax = int(0.5*(self.lmax + 1)*(self.lmax + 2)) - 1

        # Fill the rest of the last band with zeros
        temp = np.zeros(self.jmax + 1)
        temp[:len(coeffs)] = coeffs
        self.coeffs = temp

    def __add__(self, other):
        return self.coeffs + other.coeffs
        
    def __mul__(self, other):
        result = gaunt.multiply_sh_coefficients(self.coeffs, other.coeffs)
        return SHCoeffs(result)

    def plot(self, ax='shcoeffs.pdf'):
        print('test')
        if isinstance(ax, str):
            filename = ax
            f, ax0 = plt.subplots(1, 1, figsize=(4, 4))

        # Create image of spherical harmonics
        image = np.zeros((self.lmax, 2*self.lmax + 1))
        for j, c in enumerate(self.coeffs):
            l, m = util.j2lm(j)
            image[int(l/2), self.lmax + m] = c

        # Label rows and columns
        for l in range(self.lmax + 1):
            if l == 0:
                prepend = 'l='
            else:
                prepend = ''
            if l%2 == 0:
                ax0.annotate(r'$'+prepend+str(l)+'$', xy=(1, 1), xytext=(-0.75, int(l/2)),
                             textcoords='data', ha='right', va='center')

        ax0.annotate(r'$m=$', xy=(1, 1), xytext=(-0.75, -0.75),
                     textcoords='data', ha='right', va='center')
        for m in range(2*self.lmax + 1):
            ax0.annotate('$'+str(m - self.lmax)+'$', xy=(1, 1),
                         xytext=(int(m), -0.75),
                         textcoords='data', ha='center', va='center')

        # Label each pixel
        for (y,x), value in np.ndenumerate(image):
            if value != 0:
                ax0.annotate("{0:.2f}".format(value), xy=(1, 1), xytext=(x, y),
                         textcoords='data', ha='center', va='center')
            
        ax0.imshow(image, cmap='bwr', vmin=-np.max(self.coeffs), vmax=np.max(self.coeffs))
        ax0.axis('off')
        
        if isinstance(ax, str):
            f.savefig(filename)
