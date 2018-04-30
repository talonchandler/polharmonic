import numpy as np
from util import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

f, axs = plt.subplots(4, 4, figsize=(10,14),
                      gridspec_kw={'width_ratios':[1,1,1,0.05], 'hspace':0.1, 'wspace':0.25})

polar = 0
w = 2.0
out_path = 'compare-otf/'
n_px = 2**9 + 1
pad = 2**10
[X, Y] = np.meshgrid(np.linspace(-w, w, n_px),
                     np.linspace(-w, w, n_px))
R = np.sqrt(X**2 + Y**2)
Phi = np.nan_to_num(np.arctan(Y/X))

otf_funcs = [H00, H20, H22] # XXX
otf_files = ['H00an', 'H20an', 'H22an']
psf_funcs = [h00, h20, h22] # XXX
psf_files = ['hh00an', 'hh20an', 'hh22an']
ylabels = ['Image', 'OTF Profile', 'PSF Profile', 'PSF Diff']
titles = [r"${H}_0^{0(p)}(\boldsymbol{\nu})$", r"${H}_2^{0(p)}(\boldsymbol{\nu})$", r"${H}_2^{2(p)}(\boldsymbol{\nu})$"]

clabels = [[(n_px/2,0), (n_px/2,100), (n_px/2,150)],
           [(n_px/2,0), (3*n_px/4, n_px/3), (n_px/4, n_px/3), (3*n_px/4, n_px/4), (n_px/4, n_px/4)],
           [(n_px/2,0), (n_px/2,100), (n_px/2,450)]]
for i, (otf_func, otf_file) in enumerate(zip(otf_funcs, otf_files)):
    # Calculations
    otf_arr = (otf_func(R, phi=Phi, phi_p=polar)).astype('float32')
    psf_arr = (psf_funcs[i](R, phi=Phi, phi_p=polar)).astype('float32')
    xF, psf_arr_ft = myfft(otf_arr, pad=pad)

    if i==0:
        ft_norm = np.max(psf_arr_ft)
    # elif i==2: # XXX
    #     ft_norm = np.max(psf_arr_ft)

    # Save images
    # save_tiff(otf_arr, out_path+otf_file+'.tiff')
    # save_tiff(psf_arr, out_path+psf_files[i]+'.tiff')
    # save_tiff(psf_arr_ft, out_path+psf_files[i]+'ft.tiff')

    # Plot 2D OTF
    im = axs[0,i].imshow(otf_arr, cmap="bwr", vmin=-1, vmax=1)
    levels = [-0.2, -0.1, -1e-5, 1e-5, 0.1, 0.2]
    ct = axs[0,i].contour(otf_arr, levels, colors='k',linewidths=0.5)
    # plt.clabel(ct, levels, inline=1, inline_spacing=15, fmt='%1.1f', fontsize=8, manual=clabels[i])

    # Plot OTF
    axs[0,i].set_axis_off()
    x = np.linspace(-w, w, n_px)
    xd = np.linspace(-w*np.sqrt(2), w*np.sqrt(2), n_px)
    axs[1,i].plot(x, cs(otf_arr), '--', c='k', lw=1)
    axs[1,i].plot(x, cs(otf_arr, row=True), '-', c='k', lw=1)
    axs[1,i].plot(xd, ds(otf_arr), ':', c='k', lw=1)
    # axs[1,i].plot(x, cs(H20(R, phi=Phi, phi_p=None)), '-', c='g', lw=1)
    
    axs[1,i].set_xlim([-2.1,2.1])
    axs[1,i].set_ylim([-1,1])
    axs[1,i].set_xlabel(r"$\nu\, \lambda/\textrm{NA}$")
    axs[1,i].get_yaxis().set_ticks([-1, -0.5, 0.0, 0.5, 1])

    # Plot PSF
    axs[2,i].plot(x, np.abs(cs(psf_arr)), '--', c='k', lw=1)
    axs[2,i].plot(x, np.abs(cs(psf_arr, row=True)), '-', c='k', lw=1)
    # axs[2,i].plot(xd, ds(psf_arr), ':', c='k', lw=1)

    axs[2,i].plot(xF, cs(psf_arr_ft)/ft_norm, '.', c='r', lw=0.5, markersize=1)
    axs[2,i].plot(xF, cs(psf_arr_ft, row=True)/ft_norm, '.', c='r', lw=0.5, markersize=1)
    
    axs[2,i].set_xlim([-2.1,2.1])
    axs[2,i].set_ylim([-1,1])
    axs[2,i].set_xlabel(r"$r'_o \textrm{NA}/\lambda$")
    axs[2,i].get_yaxis().set_ticks([-1, -0.5, 0.0, 0.5, 1])

    # Calculate diff
    [XF, YF] = np.meshgrid(xF, xF)
    RF = np.sqrt(XF**2 + YF**2)
    PhiF = np.nan_to_num(np.arctan(YF/XF))
    psf_arrF = (psf_funcs[i](RF, phi=PhiF, phi_p=polar)).astype('float32')
        
    axs[3,i].plot(xF, cs(psf_arr_ft)/ft_norm - np.abs(cs(psf_arrF)), '--', c='k', lw=0.5)
    axs[3,i].plot(xF, cs(psf_arr_ft, row=True)/ft_norm - np.abs(cs(psf_arrF, row=True)), '-', c='k', lw=0.5)

    axs[3,i].set_xlim([-2.1,2.1])
    axs[3,i].set_ylim([-0.1,0.1])
    axs[3,i].set_xlabel(r"$r'_o \textrm{NA}/\lambda$")
    axs[3,i].get_yaxis().set_ticks([-0.1, -0.05, 0.0, 0.05, 0.1])

    axs[i,0].annotate(ylabels[i], xy=(0,0), xytext=(-0.2, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
    axs[0,i].set_title(titles[i])

# Final plotting settings
cbar = f.colorbar(im, cax=axs[0,3], orientation='vertical', ticks=[-1,-0.5, 0, 0.5, 1])
axs[0,3].get_xaxis().set_ticks([])
axs[0,3].yaxis.tick_right()
axs[1,3].set_axis_off()
axs[2,3].set_axis_off()
axs[3,3].set_axis_off()

# Save
plt.savefig('otf-compare.pdf', bbox_inches='tight')
