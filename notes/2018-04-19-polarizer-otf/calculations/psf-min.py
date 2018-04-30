import numpy as np
from util import *
import matplotlib.pyplot as plt
import matplotlib

f, axs = plt.subplots(1, 2, figsize=(6,4),
                      gridspec_kw={'width_ratios':[1,1], 'hspace':0.1, 'wspace':0.1})
w = 1.5
h = 4.5
NA = 1.0

x = np.linspace(-w, w, 2**10)

pols = [None, 0]
s = 0.65*NA

for i, pol in enumerate(pols):
    print(pol)
    center = h00(x, phi=0, NA=NA, n=1.33, phi_p=pol)
    shiftr = h00(x-s, phi=0, NA=NA, n=1.33, phi_p=pol)
    shiftl = h00(x+s, phi=0, NA=NA, n=1.33, phi_p=pol)

    axs[i].plot(x, shiftl + 3.3, '-k', lw=0.5)
    axs[i].plot(x, center + 2.3, '-k', lw=0.5)
    axs[i].plot(x, shiftr + 1.3, '-k', lw=0.5)
    axs[i].plot(x, center + shiftl + shiftr, '-k', lw=0.5)
    axs[i].plot([s,s], [-100, 100], ':k', lw=0.5)
    axs[i].plot([-s,-s], [-100, 100], ':k', lw=0.5)    
    axs[i].plot([0,0], [-100, 100], ':k', lw=0.5)

    axs[i].set_xlim([-w,w])
    axs[i].set_ylim([0,h])

axs[0].annotate(r"${h'}_0^{0(p)}(x + x')$", xy=(0,0), xytext=(-1.6, 3.3), textcoords='data', ha='right', va='center')
axs[0].annotate(r"${h'}_0^{0(p)}(x)$", xy=(0,0), xytext=(-1.6, 2.3), textcoords='data', ha='right', va='center')
axs[0].annotate(r"${h'}_0^{0(p)}(x - x')$", xy=(0,0), xytext=(-1.6, 1.3), textcoords='data', ha='right', va='center')
axs[0].annotate(r"Sum", xy=(0,0), xytext=(-1.6, 0), textcoords='data', ha='right', va='center')
axs[0].set_title(r"Without polarizer")
axs[1].set_title(r"With polarizer")

axs[0].set_axis_off()
axs[1].set_axis_off()

# Plot PSFs
plt.savefig('psf-min.pdf', bbox_inches='tight')
