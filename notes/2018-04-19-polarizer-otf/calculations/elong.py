import numpy as np
from util import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from scipy import optimize

def FWHM(X,Y):
    half_max = max(Y) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    #plot(X,d) #if you are interested
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0]
    right_idx = np.where(d < 0)[-1]
    return X[right_idx] - X[left_idx]

NAs = np.linspace(0, 1.0, 400)

def ellipse_ratio(NA):
    x = np.linspace(-2, 2, 2**12)
    ay = h00(x, phi=0, NA=NA, n=1.33, phi_p=0)
    by = h00(x, phi=np.pi/2, NA=NA, n=1.33, phi_p=0)
    a = FWHM(x, ay)
    b = FWHM(x, by)
    return a/b

ellipse_ratio_v = np.vectorize(ellipse_ratio)
ratios = ellipse_ratio_v(NAs)

f, ax = plt.subplots(1, 1, figsize=(4,4))
ax.plot(NAs, ratios, '-k')
ax.set_xlim([0, 1.0])
ax.set_ylim([1.0, 1.5])
ax.set_xlabel('NA')
ax.set_ylabel(r"${h'}_0^0$ semi-major axis/semi-minor axis")
plt.savefig('elong.pdf', bbox_inches='tight')

# polar = -np.pi/3
# w = 2.15
# out_path = 'otfs/'
# n_px = 2**10 + 1
# [X, Y] = np.meshgrid(np.linspace(-w, w, n_px),
#                      np.linspace(-w, w, n_px))
# R = np.sqrt(X**2 + Y**2)
# Phi = np.nan_to_num(np.arctan(Y/X))

# psf_funcs = [H00, H20, H22, H2_2]
# psf_files = ['H00an', 'H20an', 'H22an', 'H2_2an']
# clabels = [[(n_px/2,0), (n_px/2, n_px/10), (n_px/2,n_px/7)],
#            [(n_px/2,0), (n_px/4 + n_px/16, n_px/4), (n_px/4 + n_px/16, n_px/2 - n_px/10)],
#            [(n_px/2,0), (n_px/2, n_px/10), (n_px/2, n_px/5)],
#            [(n_px/2,0), (n_px/2, n_px/10), (n_px/2, n_px/5)]]
# for i, (psf_func, psf_file) in enumerate(zip(psf_funcs, psf_files)):
#     psf_arr = (psf_func(R, phi=Phi, phi_p=polar)).astype('float32')
#     xF, ftraw = myfft(psf_arr)
#     save_tiff(psf_arr, out_path+psf_file+'.tiff')
#     save_tiff(ftraw, out_path+psf_file+'ft.tiff')
#     im = axs[0,i].imshow(psf_arr, cmap="bwr", vmin=-1, vmax=1)
#     levels = [-0.2, -0.1, -1e-5, 1e-5, 0.1, 0.2]
#     ct = axs[0,i].contour(psf_arr, levels, colors='k',linewidths=0.5)
#     # if i!=3:
#     #     plt.clabel(ct, levels, inline=1, inline_spacing=15, fmt='%1.1f', fontsize=8, manual=clabels[i])
           
#     axs[0,i].set_axis_off()
#     x = np.linspace(-w, w, n_px)
#     xd = np.linspace(-w*np.sqrt(2), w*np.sqrt(2), n_px)
#     axs[1,i].plot(x, cs(psf_arr), '--', c='k', lw=1)
#     axs[1,i].plot(x, cs(psf_arr, row=True), '-', c='k', lw=1)
#     # axs[1,i].plot(xd, ds(psf_arr), ':', c='k', lw=1)
    
#     axs[1,i].set_xlim([-2.1,2.1])
#     axs[1,i].set_ylim([-1,1])
#     axs[1,i].set_xlabel(r"$\nu\, \lambda/\textrm{NA}$")
#     axs[1,i].get_yaxis().set_ticks([-1, -0.5, 0.0, 0.5, 1])

# axs[0,0].annotate(r"Image", xy=(0,0), xytext=(-0.25, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
# axs[1,0].annotate(r"Profiles", xy=(0,0), xytext=(-0.25, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
# axs[0,0].set_title(r"${H}_0^{0(p)}(\boldsymbol{\nu})$")
# axs[0,1].set_title(r"${H}_2^{0(p)}(\boldsymbol{\nu})$")
# axs[0,2].set_title(r"${H}_2^{2(p)}(\boldsymbol{\nu})$")
# axs[0,3].set_title(r"${H}_2^{-2(p)}(\boldsymbol{\nu})$")

# axs[0,4].get_xaxis().set_ticks([])
# axs[0,4].yaxis.tick_right()
# cbar = f.colorbar(im, cax=axs[0,4], orientation='vertical', ticks=[-1,-0.5, 0, 0.5, 1])
# axs[1,4].set_axis_off()

# # Plot PSFs
# plt.savefig('otfs.pdf', bbox_inches='tight')
