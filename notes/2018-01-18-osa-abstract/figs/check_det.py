from polharmonic import det
from sympy import Symbol, lambdify
import matplotlib.pyplot as plt
import numpy as np

theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)

dd = det.Detector(NA=1.1, theta_optical_axis=0)
dd2 = det.Detector(NA=0.8, theta_optical_axis=0)

prf = dd.det_prf()
prff = lambdify(theta, prf)
prf2 = dd2.det_prf()
prff2 = lambdify(theta, prf2)

x = np.linspace(0, 2*np.pi, 100)
y = prff(x)
y2 = prff2(x)

f, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.plot(x, y, '-k')
ax.plot(x, y2, '-b')
ax.set_xlabel('theta')
ax.set_ylabel('det prf')
ax.set_xlim([0, 2*np.pi])
ax.set_ylim([0, 1])
print(np.min(y), np.max(y))
f.savefig('det.png', dpi=300)

