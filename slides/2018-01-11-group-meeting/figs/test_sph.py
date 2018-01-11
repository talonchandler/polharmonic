from mayavi import mlab
import numpy as np
from scipy.special import sph_harm
import os

# Create a sphere
r = 0.75
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:101j, 0:2*pi:101j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
mlab.clf()

for l in range(0, 6):
    for m in range(-l, l+1):
        if m == 0:
            s = sph_harm(m, l, theta, phi).real
        if m > 0:
            s = np.sqrt(2)*((-1)**m)*(sph_harm(m, l, theta, phi).real)
        if m < 0:
            s = np.sqrt(2)*((-1)**m)*(sph_harm(np.abs(m), l, theta, phi).imag)
        n = s.clip(max=0)
        p = s.clip(min=0)*(-1)
        mlab.mesh(p*x - m, p*z - 1.5*l, p*y,
                  color=(1, 0, 0), representation='surface')
        mlab.mesh(n*x - m, n*z - 1.5*l, n*y,
                  color=(0, 0, 1), representation='surface')
        
mlab.view(azimuth=0, elevation=0, distance=20, focalpoint=None,
          roll=None, reset_roll=True, figure=None)
#mlab.show()
mlab.savefig('picture.png', magnification=10)
