import numpy as np
from util import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

s = 0.75
f, axs = plt.subplots(2, 5, figsize=(s*14,s*6),
                      gridspec_kw={'width_ratios':[1,1,1,1,0.05], 'hspace':0.03, 'wspace':0.30})

polar = -np.pi/2
w = 1.25
out_path = 'psfs/'
n_px = 1024
[X, Y] = np.meshgrid(np.linspace(-w, w, n_px),
                     np.linspace(-w, w, n_px))
R = np.sqrt(X**2 + Y**2)
Phi = np.nan_to_num(np.arctan(Y/X))

psf_funcs = [h00, h20, h22, h2_2]
psf_files = ['hh00an', 'hh20an', 'h22an', 'h2_2an']
clabels = [[(n_px/2, n_px/3)], [(n_px/2,n_px/3), (n_px/4,n_px/2), (3*n_px/4,n_px/2)], [(n_px/2, n_px/3)]]
for i, (psf_func, psf_file) in enumerate(zip(psf_funcs, psf_files)):
    psf_arr = (psf_func(R, phi=Phi, phi_p=polar)).astype('float32')
    #save_tiff(psf_arr, out_path+psf_file+'.tiff')
    levels = [-0.1, 0.1]
    im = axs[0,i].imshow(psf_arr, cmap="bwr", vmin=-1, vmax=1)
    ct = axs[0,i].contour(psf_arr, levels, colors='k',linewidths=0.5)
    # plt.clabel(ct, levels, inline=1, inline_spacing=15, fmt='%1.1f', fontsize=8, manual=clabels[i])
    axs[0,i].set_axis_off()
    x = np.linspace(-w, w, n_px)
    axs[1,i].plot(x, cs(psf_arr), '--', c='k', lw=1)
    axs[1,i].plot(x, cs(psf_arr, row=True), '-', c='k', lw=1)
    axs[1,i].set_xlim([-w,w])
    axs[1,i].set_ylim([-1,1])
    axs[1,i].set_xlabel(r"$r_o'\, \textrm{NA}/\lambda$")
    axs[1,i].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])

axs[0,0].annotate(r"Image", xy=(0,0), xytext=(-0.25, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
axs[1,0].annotate(r"Profiles", xy=(0,0), xytext=(-0.25, 0.5), textcoords='axes fraction', ha='center', va='center', rotation=90)
axs[0,0].set_title(r"${h'}_0^{0(p)}(\mathbf{r}_o)$")
axs[0,1].set_title(r"${h'}_2^{0(p)}(\mathbf{r}_o)$")
axs[0,2].set_title(r"${h'}_2^{2(p)}(\mathbf{r}_o)$")
axs[0,3].set_title(r"${h'}_2^{-2(p)}(\mathbf{r}_o)$")

axs[0,4].get_xaxis().set_ticks([])
axs[0,4].yaxis.tick_right()
cbar = f.colorbar(im, cax=axs[0,4], orientation='vertical', ticks=[-1,-0.5, 0, 0.5, 1])
axs[1,4].set_axis_off()

# Plot PSFs
plt.savefig('psfs.pdf', bbox_inches='tight')
