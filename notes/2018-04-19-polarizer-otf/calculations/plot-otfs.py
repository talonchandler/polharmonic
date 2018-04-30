import numpy as np
from util import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

s = 0.75
f, axs = plt.subplots(2, 5, figsize=(s*14,s*6),
                      gridspec_kw={'width_ratios':[1,1,1,1,0.05], 'hspace':0.03, 'wspace':0.3})

polar = -np.pi/3
w = 2.15
out_path = 'otfs/'
n_px = 2**10 + 1
[X, Y] = np.meshgrid(np.linspace(-w, w, n_px),
                     np.linspace(-w, w, n_px))
R = np.sqrt(X**2 + Y**2)
Phi = np.nan_to_num(np.arctan(Y/X))

psf_funcs = [H00, H20, H22, H2_2]
psf_files = ['H00an', 'H20an', 'H22an', 'H2_2an']
clabels = [[(n_px/2,0), (n_px/2, n_px/10), (n_px/2,n_px/7)],
           [(n_px/2,0), (n_px/4 + n_px/16, n_px/4), (n_px/4 + n_px/16, n_px/2 - n_px/10)],
           [(n_px/2,0), (n_px/2, n_px/10), (n_px/2, n_px/5)],
           [(n_px/2,0), (n_px/2, n_px/10), (n_px/2, n_px/5)]]
for i, (psf_func, psf_file) in enumerate(zip(psf_funcs, psf_files)):
    psf_arr = (psf_func(R, phi=Phi, phi_p=polar)).astype('float32')
    xF, ftraw = myfft(psf_arr)
    save_tiff(psf_arr, out_path+psf_file+'.tiff')
    save_tiff(ftraw, out_path+psf_file+'ft.tiff')
    im = axs[0,i].imshow(psf_arr, cmap="bwr", vmin=-1, vmax=1)
    levels = [-0.2, -0.1, -1e-5, 1e-5, 0.1, 0.2]
    ct = axs[0,i].contour(psf_arr, levels, colors='k',linewidths=0.5)
    # if i!=3:
    #     plt.clabel(ct, levels, inline=1, inline_spacing=15, fmt='%1.1f', fontsize=8, manual=clabels[i])
           
    axs[0,i].set_axis_off()
    x = np.linspace(-w, w, n_px)
    xd = np.linspace(-w*np.sqrt(2), w*np.sqrt(2), n_px)
    axs[1,i].plot(x, cs(psf_arr), '--', c='k', lw=1)
    axs[1,i].plot(x, cs(psf_arr, row=True), '-', c='k', lw=1)
    # axs[1,i].plot(xd, ds(psf_arr), ':', c='k', lw=1)
    
    axs[1,i].set_xlim([-2.1,2.1])
    axs[1,i].set_ylim([-1,1])
    axs[1,i].set_xlabel(r"$\nu\, \lambda/\textrm{NA}$")
    axs[1,i].get_yaxis().set_ticks([-1, -0.5, 0.0, 0.5, 1])

axs[0,0].annotate(r"Image", xy=(0,0), xytext=(-0.25, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
axs[1,0].annotate(r"Profiles", xy=(0,0), xytext=(-0.25, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
axs[0,0].set_title(r"${H}_0^{0(p)}(\boldsymbol{\nu})$")
axs[0,1].set_title(r"${H}_2^{0(p)}(\boldsymbol{\nu})$")
axs[0,2].set_title(r"${H}_2^{2(p)}(\boldsymbol{\nu})$")
axs[0,3].set_title(r"${H}_2^{-2(p)}(\boldsymbol{\nu})$")

axs[0,4].get_xaxis().set_ticks([])
axs[0,4].yaxis.tick_right()
cbar = f.colorbar(im, cax=axs[0,4], orientation='vertical', ticks=[-1,-0.5, 0, 0.5, 1])
axs[1,4].set_axis_off()

# Plot PSFs
plt.savefig('otfs.pdf', bbox_inches='tight')
