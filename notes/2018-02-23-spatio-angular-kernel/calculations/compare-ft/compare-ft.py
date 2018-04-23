import numpy as np
from util import *
import matplotlib.pyplot as plt

# Setup
out = 'output/'
NA = 0.8
n = 1.33
n_pixels = 257
[X, Y] = np.meshgrid(np.linspace(-2, 2, n_pixels),
                     np.linspace(-2, 2, n_pixels))
R = np.sqrt(X**2 + Y**2)

# Analytic OTF
Hxan = (A(R)).astype('float32')
Hzan = (B(R, NA=NA, n=n)).astype('float32')
H00an = (H00(R, NA=NA, n=n)).astype('float32')
H20an = (H20(R, NA=NA, n=n)).astype('float32')
save_tiff(Hxan, out+'Hxan.tiff')
save_tiff(Hzan, out+'Hzan.tiff')
save_tiff(H00an, out+'H00an.tiff')
save_tiff(H20an, out+'H20an.tiff')

# Analytic PSF
hhxan = (a(R)).astype('float32')
hhzan = (b(R, NA=NA, n=n)).astype('float32')
hh00an = (h00(R, NA=NA, n=n)).astype('float32')
hh20an = (h20(R, NA=NA, n=n)).astype('float32')
save_tiff(hhxan, out+'hhxan.tiff')
save_tiff(hhzan, out+'hhzan.tiff')
save_tiff(hh00an, out+'hh00an.tiff')
save_tiff(hh20an, out+'hh20an.tiff')

# Computed PSF
xF, hhxftraw = myfft(Hxan)
xF, hhzftraw = myfft(Hzan)
xF, hh00ftraw = myfft(H00an)
xF, hh20ftraw = myfft(H20an)

hhxft = hhxftraw/np.max(hhxftraw)
hhzft = hhzftraw/np.max(hhxftraw)
hh00ft = hh00ftraw/np.max(hh00ftraw)
hh20ft = hh20ftraw/np.max(hh00ftraw)

save_tiff(hhxft, out+'hhxft.tiff')
save_tiff(hhzft, out+'hhzft.tiff')
save_tiff(hh00ft, out+'hh00ft.tiff')
save_tiff(hh20ft, out+'hh20ft.tiff')

# Plotting

# Plot PSFs
x = np.linspace(-2, 2, n_pixels)
f, axs = plt.subplots(1, 2, figsize=(10,5))

axs[0].plot(x, cs(hh00an), '-', c='k', lw=1, label=r"${h_0^0}^{(p)}(r_o')$")
#axs[0].plot(xF, cs(hh00ft), ':', c='k',lw=1, label=r"$h_0^0(r_d')$ FT")

axs[0].plot(x, cs(hh20an), '-', c='b',lw=1, label=r"${h_2^0}^{(p)}(r_o')$")
#axs[0].plot(xF, cs(hh20ft), ':', c='b',lw=1, label=r"$h_0^2(r_d')$ FT")

axs[0].plot(x, cs(hhxan), '-', c='r',lw=1, label=r"$h_x^{(p)}(r_o')$")
#axs[0].plot(xF, cs(hhxft), ':', c='r',lw=1, label=r"$h_x(r_d')$ FT")

axs[0].plot(x, cs(hhzan), '-', c='g',lw=1, label=r"$h_z^{(p)}(r_o')$")
#axs[0].plot(xF, cs(hhzft), ':', c='g',lw=1, label=r"$h_z(r_d')$ FT")

# Ax 0 setup
axs[0].set_xlim([-2, 2])
axs[0].set_ylim([-1, 1])
axs[0].set_xlabel(r"$r_o'\, \textrm{NA}/\lambda$")
axs[0].set_title(r'')
axs[0].get_xaxis().set_ticks([-2, -1, 0, 1, 2])
axs[0].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
axs[0].xaxis.labelpad = 10
axs[0].legend(frameon=False, loc=4)

# Plot OTFS
x = np.linspace(-2, 2, n_pixels)


axs[1].plot(x, cs(H00an), '-', c='k', lw=1, label=r"${H_0^0}^{(p)}(\nu)$")
axs[1].plot(x, cs(H20an), '-', c='b',lw=1, label=r"${H_2^0}^{(p)}(\nu)$")
axs[1].plot(x, cs(Hxan), '-', c='r',lw=1, label=r"$H_x^{(p)}(\nu)$")
axs[1].plot(x, cs(Hzan), '-', c='g',lw=1, label=r"$H_z^{(p)}(\nu)$")
N = (1/np.pi)*((NA/n)**2)
poly = (1.0/3.0)*(13.0 - 10.0*(np.abs(x/2)**2))

axs[1].set_xlim([-1, 1])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel(r"$\nu\lambda/\textrm{NA}$")
axs[1].set_title(r'')
axs[1].get_xaxis().set_ticks([-2, -1, 0, 1, 2])
axs[1].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
axs[1].xaxis.labelpad = 10
axs[1].legend(frameon=False, loc=4)

plt.savefig('psfs.pdf')

