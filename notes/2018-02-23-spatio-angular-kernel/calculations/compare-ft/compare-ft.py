import numpy as np
from util import *
import matplotlib.pyplot as plt

# Setup
out = 'output/'
NA = 0.8
n = 1.33
n_pixels = 256
[X, Y] = np.meshgrid(np.linspace(-2, 2, n_pixels),
                     np.linspace(-2, 2, n_pixels))
R = np.sqrt(X**2 + Y**2)

# Analytic OTF
Han = (A(R)).astype('float32')
H00an = (H00(R)).astype('float32')
H20an = (H20(R)).astype('float32')
save_tiff(Han, out+'Han.tiff')
save_tiff(H00an, out+'H00an.tiff')
save_tiff(H20an, out+'H20an.tiff')

# Analytic PSF
hhan = (a(R)).astype('float32')
hh00an = (h00(R)).astype('float32')
hh20an = (h20(R)).astype('float32')
save_tiff(hhan, out+'hhan.tiff')
save_tiff(hh00an, out+'hh00an.tiff')
save_tiff(hh20an, out+'hh20an.tiff')

# Computed PSF
xF, hhft = myfft(Han)
xF, hh00ft = myfft(H00an)
xF, hh20ft = myfft(H20an)

save_tiff(hhft, out+'hhft.tiff')
save_tiff(hh00ft, out+'hh00ft.tiff')
save_tiff(hh20ft, out+'hh20ft.tiff')

# Plotting

# Plot PSFs
x = np.linspace(-2, 2, n_pixels)
y2 = h20(x, NA=NA)
y3 = h20(x, NA=NA)*(-np.sqrt(5))
f, axs = plt.subplots(1, 2, figsize=(10,5))

axs[0].plot(x, cs(hh00an), '-', c='k', lw=1, label=r"$h_0^0(r_d')$")
axs[0].plot(xF, cs(hh00ft)/np.max(hh00ft), ':', c='k',lw=1, label=r"$h_0^0(r_d')$ FT")

axs[0].plot(x, cs(hh20an), '-', c='b',lw=1, label=r"$h_0^2(r_d')$")
axs[0].plot(xF, cs(hh20ft)/np.max(hh00ft), ':', c='b',lw=1, label=r"$h_0^2(r_d')$ FT")

# axs[0].plot(x, y3, ':', c='k',lw=1, label=r"Normalized $h_0^2(\nu)$")
# axs[0].plot(x, H00ft/np.max(H00ft), ':', c='b',lw=1, label=r"H00ft calc")
# axs[0].plot(x, H20ft/np.max(H00ft), ':', c='r',lw=1, label=r"H20ft calc")

axs[0].plot(x, a(x), '-', c='r',lw=1, label=r"Airy")
axs[0].plot(xF, cs(hhft)/np.max(hhft), ':', c='r',lw=1, label=r"Airy FT")

axs[0].set_xlim([-2, 2])
axs[0].set_ylim([-1, 1])
axs[0].set_xlabel(r"$r_d'\, \textrm{NA}/(M\lambda)$")
axs[0].set_title(r'Spatio-angular PSF')
axs[0].get_xaxis().set_ticks([-2, -1, 0, 1, 2])
axs[0].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
axs[0].xaxis.labelpad = 10
axs[0].legend(frameon=False, loc=8)

# Plot OTFS
x = np.linspace(-2, 2, n_pixels)
y2 = H20(x, NA=NA, n=n)#*(np.sqrt(5)*(1 + 0.5*(NA/n)**2))/(-1 + ((NA/n)**2))
y3 = H20(x, NA=NA, n=n)*(np.sqrt(5)*(1 + 0.5*(NA/n)**2))/(-1 + ((NA/n)**2))
y4 = A(x)

axs[1].plot(x, cs(H00an), '-', c='k', lw=1, label=r"$H_0^0(r_d')$")
axs[1].plot(x, cs(H20an), '-', c='b',lw=1, label=r"$H_0^2(\nu)$")
axs[1].plot(x, cs(Han), '-', c='r',lw=1, label=r"Airy OTF")

# axs[1].plot(x, y3, ':', c='k',lw=1, label=r"Normalized $H_0^2(\nu)$")
# axs[1].plot(x, H00im, ':', c='b',lw=1, label=r"H00 calc")
# axs[1].plot(x, H20im, ':', c='r',lw=1, label=r"H20 calc")

axs[1].set_xlim([-1, 1])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel(r"$\nu\, (M\lambda)/\textrm{NA}$")
axs[1].set_title(r'Spatio-angular OTF')
axs[1].get_xaxis().set_ticks([-2, -1, 0, 1, 2])
axs[1].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
axs[1].xaxis.labelpad = 10
axs[1].legend(frameon=False, loc=8)

plt.savefig('psfs.pdf')

